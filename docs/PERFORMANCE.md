# Performance vs SciPy

Goal: core hot paths should be **as fast as SciPy where the comparison is fair**,
and **faster where signal4s architecture allows** (plans, workspaces, sparse
structure), without regressing fixture parity or laws.

## Fairness rules

1. **Same problem.** Identical `N`, `M`, section count, `up`/`down`, and numeric
   inputs (or the same deterministic generator).
2. **Cold and amortized paths.** Report one-shot FFT convolution against
   `signal.convolve(..., method="fft")`, and compare reused `ConvolutionPlan`
   against an explicit SciPy `rfft(x) * precomputedKernelSpectrum` peer. Reused
   `FftPlan` + workspace and `SosRunner.processInto` likewise need repeated-call
   peers. FFT and OLA plans reuse mutable work buffers and are single-owner (not
   concurrent), matching filter runners.
3. **Explicit methods.** Time `Direct` and `Fft` separately; do not hide Auto
   dispatch inside a vs-SciPy cell.
4. **Warm-up.** Discard JIT / NumPy warm-up iterations before the timed window.
   Report the **median of several trials**.
5. **Correctness gate.** Existing SciPy fixtures, FFT↔direct parity across a
   size grid, and law suites must stay green. Speed never overrides numeric
   contracts in [`SEMANTICS.md`](SEMANTICS.md).
6. **No narrow benchmaxxing.** Optimize algorithms that help many sizes (e.g.
   5-smooth FFT lengths, packed real FFT when power-of-two). Do not special-case
   the exact receipt `(N,M)` pairs or game the harness.

## Targets (JVM, after warm-up)

| Op | Fair SciPy peer | Target ratio (signal4s / SciPy) | Stretch |
|---|---|---|---|
| Direct convolve full | `signal.convolve(..., method="direct")` | ≤ 1.5× | ≤ 1.0× |
| Direct convolve valid | `signal.convolve(..., mode="valid", method="direct")` | ≤ 1.5× | ≤ 1.0× |
| Causal FIR / Input(Zero) | `signal.lfilter(h, [1], x)` | ≤ 1.5× | ≤ 1.0× |
| FFT convolve (plan reuse) | precomputed-kernel `rfft`/multiply/`irfft` | ≤ 1.5× | ≤ 1.0× |
| Circular FFT convolve | `rfft`/multiply/`irfft` of a periodized kernel | ≤ 1.5× | ≤ 1.0× |
| Power-of-two FFT (workspace) | `numpy.fft.fft` | ≤ 2.0× | ≤ 1.2× |
| SOS filter (runner reuse) | `signal.sosfilt` | ≤ 1.5× | ≤ 1.0× |
| upfirdn | `signal.upfirdn` | ≤ 1.5× | ≤ 1.0× |

Ratios > 1 mean signal4s is slower. Report both `ms/call` and the ratio in
[`../benchmarks/receipts/vs-scipy.md`](../benchmarks/receipts/vs-scipy.md).

Raw complex power-of-two and 2/3/5-smooth transforms, plus even-length real
transforms, use the JTransforms JVM backend. The real backend uses
JTransforms' packed layout directly and avoids routing real input through a
split-complex transform.
On macOS, dense direct FIRs use the system `vDSP_convD` routine when its
Accelerate symbol is available. Other operating systems retain the portable
scalar and Vector API kernels. The macOS backend requires JDK 22 or later and
the configured `--enable-native-access=ALL-UNNAMED` runtime flag.

FFT convolution selects the next **power of two** (`FastFftLength`) so the packed
real FFT path applies. The portable mixed-radix (2/3/5-smooth) engine is
correct, uses pooled part buffers and specialized radix-2/4 combines, and is
covered by a DFT grid. It remains the Scala.js/fallback path and is still
roughly 2–6× slower per transform than a JVM packed-real transform at the next
power of two, so it is not selected for convolution lengths yet.

## Harness

```bash
# SciPy side (pinned env from fixtures/generate/requirements.txt)
python3 fixtures/generate/bench_vs_scipy.py

# signal4s side (JVM)
sbt "fftJVM / Test / runMain signal4s.fft.bench.VsScipyBench"
```

Both write machine-readable JSON under `benchmarks/receipts/`. The Scala main
merges them into `vs-scipy.md`.

## Optimization discipline

1. Keep a scalar / sparse-oracle path for laws when specializing kernels.
2. Prefer specialized `while` loops and precomputed tables over framework layers.
3. Measure before and after on the vs-SciPy shapes; update the receipt.
4. Do not change public semantics to chase a ratio.

## Priority order

1. Comparative harness + baseline receipt *(done)*
2. Direct convolution Zero / causal specialization *(done; i-blocked scatter,
   x86 Vector API, and macOS vDSP)*
3. upfirdn without materializing the full upsample buffer *(done; decimate tap
   banks and a division-free general polyphase recurrence)*
4. FFT twiddle tables, packed real FFT, Hermitian multiply *(done)*
5. JVM JTransforms backend for power-of-two and smooth complex FFT *(done; closes FFT-conv ≤1.5×)*
6. Verified Stockham autosort *(correct; slower than bit-rev here — tests only)*
7. JVM Vector API FIR on x86 AVX *(landed; disabled on aarch64 where width=2 / emulated)*
8. macOS vDSP full FIR *(done; one `vDSP_convD` call, off-heap reusable staging)*
9. Auto cost-model validation against measured costs

### Hot-path profile (A1/B1)

- **Direct:** `scatterFull` is ~95–98% of e2e at large `N`; adopt ~1%. On
  macOS, dense full FIRs (`N ≥ 64`, `16 ≤ M ≤ 2044`) use one `vDSP_convD`
  call with thread-local off-heap staging. The vDSP-eligible receipt cells are
  0.18–0.48× their SciPy peers on Apple Silicon. The current receipt table also
  contains a 0.57× Direct cell (`N=8192`, `M=128`), so these are scoped sample
  results rather than a universal bound. Other platforms use the scalar/
  i-blocked or x86 Vector API kernels.
- **Causal FIR:** `Input(Zero)` with a causal kernel now reuses the native/vector
  full-scatter tier when \(M \le N\). This preserves scalar work when a long
  kernel would make the full tail wasteful; the measured 1024–8192 grid is
  0.21–0.48× SciPy `lfilter`.
- **Valid direct convolution:** native/vector backends similarly form one full
  convolution and retain the fully-overlapped interval. This changes the
  1024–8192 causal-kernel grid from 1.76–4.09× to 0.36–0.51× SciPy direct
  `mode="valid"` while retaining the scalar path where it is more appropriate.
- **Circular convolution:** when `signal4s-fft` is installed, the core API
  automatically selects a packed-real FFT for native-backed even periods and
  retains direct work on the portable conservative crossover. The current
  1024–8192 circular FFT grid is 0.62–1.03× the equivalent NumPy transform.
- **Smooth complex FFT:** the JVM now routes 2/3/5-smooth non-power-of-two
  transforms through JTransforms while retaining portable mixed-radix as a
  tested fallback. On the workspace bench, 2160 improved from 0.157 to 0.069
  ms/call and 8640 from 0.440 to 0.091 ms/call.
- **upfirdn:** the general `up > 1`, `down > 1` path tracks the quotient and
  remainder of `inputIndex * up / down`. Its inner loop increments contiguous
  output positions instead of dividing once per retained tap. The receipt grid
  is 0.48–1.42× SciPy in the current receipt grid without materializing the
  upsampled signal. The range is a measured snapshot, not a platform guarantee.
- **FFT-conv:** with the packed JTransforms real backend, the symmetric planned
  grid is ≤1.0× SciPy (8192/128 is ~0.96×). One-shot timings are reported
  separately because both implementations rebuild the kernel transform. Planned
  `Valid` and `Input(Zero)` calls copy only the requested post-IFFT region into
  a region-sized result buffer.
- **Real FFT:** workspace-reused `rfft` is 0.07–1.25× NumPy over lengths
  64–8192. Even-length `RealFftPlan`, Welch, and STFT share this backend.
- **Overlap-add:** plans retain FFT and block buffers plus a region-sized
  accumulation buffer. Each apply allocates only the immutable requested result,
  not a full intermediate result or a region slice.

Receipts use a full settle pass, then the **median of 7 trials** (100 iters,
50 warm-up) on both sides. They are machine-, JDK-, and platform-specific;
the current receipt is dated 2026-08-01 and E11 performance work remains
active.
