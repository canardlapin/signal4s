# signal4s — Architecture Proposal

**Status:** Proceed
**Working name:** `signal4s`
**Date captured:** 2026-07-31

---

## Verdict

Build a compact calculus of sampled signals and systems. Do **not** build a Scala translation of `scipy.signal` or `gsignal`.

`gsignal` is a coverage checklist, not an architectural template. SciPy is a numerical reference, not an ontology. Both treat a signal as an array and select semantics with strings. `signal4s` should win on precision of concepts, explicit execution, and laws that make extensions trustworthy—not on breadth.

### Conceptual spine

> **A signal has coordinates. A system has semantics. A representation describes the system. A plan selects an algorithm. A runner owns state.**

### Load-bearing architecture

```text
Signal + Kernel + System Representation + Plan + Runner + Typed Result Axes
```

### First slice (highest leverage)

1. `Kernel` with explicit zero-lag origin
2. Semantically exact convolution regions and boundaries
3. Direct convolution
4. Gale `ConvolutionOperator` with a proven adjoint
5. Laws that make future implementations interchangeable

---

## 1. Idea (one paragraph)

A discrete signal is a function \(\mathbb Z \to \mathbb R\) of which a computer holds a finite observation. LTI systems are polynomials (FIR) or rational functions (IIR) of the delay operator. Convolution, filtering, Fourier analysis, spectral estimation, and resampling all sit on five organizing pieces: finite support, delay/LTI systems, state realization, transform plans, and lattice changes. The library makes axes, representations, algorithms, and mutable state into distinct values, so common operations stay inspectable and interchangeable under floating-point laws.

---

## 2. Design constitution (8 rules)

1. **Axes are data.** Sample rate, origin, lags, frame times, frequency bins, and latency live on values, not only in docs.
2. **Representations are distinct values.** TF, ZPK, SOS, and state-space are not interchangeable matrix shapes.
3. **Meaning is separate from method.** Direct, FFT, and overlap-add return the same semantic result; the algorithm is a plan.
4. **State is explicit.** Filter descriptions are immutable; streaming uses a single-owner runner.
5. **Allocation has a visible boundary.** Allocating conveniences for ordinary use; destination-writing / workspace reuse for hot paths.
6. **Physical conventions are never implicit.** No unexplained `Double` that might mean Hz, cycles/sample, Nyquist units, or rad/sample.
7. **Mathematics supplies laws, not ceremony.** Users do not need category theory; laws govern correctness of new implementations.
8. **Start with `Double`.** Add `Float` only when a measured workload justifies it.

---

## 3. Dependencies and module boundary

| Need | Representation |
|---|---|
| Samples and taps | Gale `DVec` |
| Mutable destinations | Gale `MutableDVec` |
| Immutable construction | Gale `DVecBuilder` |
| State-space systems | Gale `DMat` |
| Finite batch convolution / deconvolution | Gale `DoubleLinearOperator` |
| Least-squares FIR, companion eigenvalues | Gale solvers / eig |
| N-D batches / transform-along-axis | Optional Ravel adapter |
| Complex FFT storage | Focused split-complex type in `signal4s` |

- **Gale:** real dependency of `signal4s-core`.
- **Ravel:** optional `signal4s-ravel` only; do not put N-D ontology in the 1-D core (Ravel is pre-1.0).
- **Gale enhancement (deferred):** if accessor overhead matters (esp. Scala.js), add a narrow strided-access SPI—do not expose backing arrays, and do not add this preemptively.

---

## 4. Itemized architecture

### 4.1 Central value types

- [ ] `SampleRate`, `Frequency`, `Seconds`, `Sampling`
- [ ] `Signal(samples: DVec, sampling: Sampling)`
- [ ] `Kernel(taps: DVec, zeroLagIndex: Int)` with `causal` / `at` / `centeredOdd` constructors
- [ ] Never silently center an even-length kernel; origin must be explicit
- [ ] Output time axis: \(t_{\mathrm{out}} = t_{\mathrm{in}} - \texttt{zeroLagIndex}/f_s\)

### 4.2 Convolution (three separate dimensions)

Separate: output support, kernel alignment, boundary, algorithm.

- [ ] `Boundary`: Zero, Constant, Clamp, Reflect, Symmetric
- [ ] `OutputRegion`: Full, Valid, Input(boundary) — call it `Input`, not `Same`
- [x] `ConvolutionMethod`: Auto, Direct, Fft, OverlapAdd(blockLength)
- [x] Circular convolution as a separate operation on \(\mathbb Z/N\mathbb Z\)
- [x] `Convolution.apply` / `Convolution.plan` → `ConvolutionPlan` (expose selected algorithm; reuse transformed kernel)
- [ ] `Correlation` returns `LagAxis` + values; normalization as enum (raw / biased / unbiased / coefficient)

Convention:

\[
r_{xy}[\ell] = \sum_n x[n]\, y[n+\ell]
\]

### 4.3 Finite convolution as Gale operator

- [ ] `Convolution.operator(...)` → `DoubleLinearOperator` with adjoint
- [ ] Law: \(\langle Ax,y\rangle \approx \langle x,A^*y\rangle\)
- [ ] Initially: zero-extension regions only; do not advertise adjoint for reflective/clamped until accumulation is correct
- [ ] Streaming runners must **not** extend `DoubleLinearOperator`

Enables: deconvolution, regularized inverse, convolutional LS, matched filtering, gradients, composition with other Gale operators.

### 4.4 Filters as systems (not coefficient arrays)

Representations:

- [ ] `Fir(kernel)`
- [ ] `DigitalTransferFunction(feedForward, feedback)` — delay powers \(z^{-k}\)
- [ ] `Biquad` with named coeffs, \(a_0 = 1\)
- [ ] `SecondOrderCascade(gain, sections)`
- [ ] `DiscreteStateSpace(a,b,c,d)`
- [ ] Analog vs digital as distinct types (or phantom `Domain`)
- [ ] SOS = normal high-order IIR realization (DF-II transposed)
- [ ] State-space for MIMO, init, composition, future Kalman work

### 4.5 Description → realization → execution

- [ ] Immutable system description; no hidden delay registers
- [ ] `filter.process(signal)` batch convenience
- [ ] `filter.newRunner()` streaming; `process` / `processInto` / `reset` / snapshot
- [ ] Document: thread-safety (generally not), latency, state shape, flush
- [ ] Chunking law: `run(x ++ y, s0) ≈ run(y, s_x)` after running `x` from `s0`
- [ ] Do **not** invent one universal `Processor[...]` hierarchy
- [ ] Zero-phase / forward–backward as separate offline algorithm (`ZeroPhase.filter`), not a runner flag

### 4.6 Frequency types

- [ ] `SampleRate`, `Frequency` (Hz), `RadiansPerSample`, `CyclesPerSample`
- [ ] Conversions require explicit sample rate; reject above Nyquist
- [x] Filter design APIs use named physical quantities, never bare `cutoff: Double`
- [ ] Design report records convention used

### 4.7 FFT and complex data

- [ ] `Complex` (scalar) + split `ComplexVector(real, imaginary)`
- [ ] `FftNormalization`: Backward / Forward / Orthonormal
- [ ] Immutable thread-safe `RealFftPlan`; mutable `RealFftWorkspace` (not cached in plan)
- [ ] Mixed-radix + Bluestein; real↔half-complex; identical JVM / Scala.js semantics
- [ ] Spectra carry axes + scaling: `RealSpectrum`, `MagnitudeSpectrum`, `PowerSpectrum`, `PowerSpectralDensity`, `CrossSpectrum`
- [ ] No bare `(frequencies, values)` tuples

### 4.8 Windows and frames

- [x] `WindowSpec` + `WindowConvention` (Symmetric vs Periodic) — always explicit
- [x] `Window` carries coherent gain, power gain, ENBW
- [x] Shared `FramePlan` (length, hop, alignment, boundary) for STFT, Welch, OLA, features

### 4.9 STFT as analysis/synthesis plan

- [x] `StftPlan` with analysis/synthesis windows, frames, FFT, sample rate
- [x] Validate invertible frame; `CanonicalDual` computed once
- [x] `TimeFrequency` with named frame/frequency axes
- [x] Spectrogram = derived view of STFT (no reimplementation)

### 4.10 Spectral estimation

**Initial:** Periodogram, Welch, CrossSpectrum
**Next:** Coherence, TransferEstimate

- [x] `WelchPlan`: frames, window, FFT, detrend, sides, scaling, average
- [x] Enums for detrending / sides / scaling / average
- [x] `WelchResult` retains segment count, effective window power, optional DoF, diagnostics

### 4.11 Multirate

- [x] Foundational op: `upfirdn` via polyphase
- [x] `RateRatio(up, down)` GCD-reduced, positive
- [x] `PolyphaseResampler` exposes rates, phase, delay, consume/produce, flush
- [x] `decimate` / `resample` compile to this substrate
- [x] Rate-changing streaming contract lives here—not bolted onto simple filter runners

### 4.12 Filter design

**FIR (initial):** windowed-sinc L/H/BP/BS; Kaiser order; maybe LS later
**IIR (initial):** Butterworth → transform → prewarped bilinear → ZPK → SOS pairing + reports

- [x] `DesignedIir(zpk, sos, report)` — SOS primary for execution
- [x] Named conversions: `Bilinear.transform`, `Sos.fromZerosPolesGain` (StateSpace deferred)
- [ ] Deferred: Chebyshev, elliptic, Bessel, Parks–McClellan, arbitrary-response, `StateSpace.fromTransferFunction`

---

## 5. Module structure

| Module | Responsibility | Cross-compile |
|---|---|---|
| `signal4s-core` | Signals, axes, kernels, direct conv/corr, FIR, TF, biquads, SOS runners, windows, errors | JVM + JS |
| `signal4s-fft` | Complex, FFT, FFT conv, framing, STFT, periodogram, Welch | JVM + JS |
| `signal4s-design` | FIR design, analog prototypes, bilinear, ZPK/SOS/TF/SS | JVM + JS |
| `signal4s-laws` | Reusable MUnit/ScalaCheck law bundles | JVM + JS |
| `signal4s-ravel` | Optional N-D / axis adapters | JVM + JS |
| `signal4s-backend-jvm-vector` | Optional Vector API kernels | later |
| `signal4s-backend-jvm-native` | Optional native FFT | later |
| `signal4s-fs2` | Optional stream adapter | later |

Layout sketch:

```text
modules/
  core/shared/src/main/scala/signal4s/
  fft/shared/src/main/scala/signal4s/fft/
  design/shared/src/main/scala/signal4s/design/
  ravel/
  laws/
  benchmarks/
```

---

## 6. Scope

### Required before 1.0

- Direct, planned FFT, overlap-add, and circular convolution
- Cross-/auto-correlation with lag axes
- FIR + TF + SOS filtering with initial/final state
- Forward–backward zero-phase filtering
- FFT / iFFT / real FFT
- Core windows (periodic vs symmetric)
- Framing; STFT + exact inverse STFT
- Periodogram and Welch PSD
- Frequency, impulse, step response; group delay
- Mean and linear detrending
- Windowed-sinc FIR design; Butterworth → SOS
- `upfirdn` and rational polyphase resampling
- Gale convolution operator with tested adjoint
- Optional Ravel axis adapters

### Explicitly deferred

Exotic windows; Chebyshev/elliptic; Parks–McClellan; wavelets; AR/ME spectral estimators; peak finding; chirps/comms generators; audio codecs/I/O; image filters; Matlab/R/SciPy API compatibility; lazy signal graphs; Cats Effect/FS2 in core; autodiff; GPU.

---

## 7. Performance model

1. **Pure surface, mutable interior** — allocating APIs + `*Into` / workspace forms.
2. **Plans amortize structure** — factorization, twiddles, transformed kernels, polyphase, window metrics.
3. **Workspaces/runners are single-owner** — never escape into immutable results.
4. **`Auto` dispatch is measured** per platform (JVM portable / Vector / native; Scala.js Node / browser); no folklore thresholds.
5. **Portable Scala FFT first**; native optional and capability-driven; must not change normalization, ordering, axes, errors, ownership, or diagnostics.
6. **Specialize for `Double`**; no iterators/`map`/`zip`/`Vector[Double]` in numerical kernels.

---

## 8. Errors and diagnostics

Sealed `SignalError` algebra (invalid rates/frequencies, Nyquist, kernel origin, length mismatch, filter coeffs/state, non-invertible frame, rate ratio, transform length, numerical failure).

Gale discipline:

- smart constructors → `Either`
- predictable numerical failure → `Either`
- qualified success → diagnostics on result
- `.orThrow` for examples/tests
- invariant-preserving ops need not wrap in `Either`

Design reports (e.g. `IirDesignReport`) record order, frequencies, pole magnitude, stability, section ordering, warnings.

---

## 9. Law suite (`signal4s-laws`)

| Area | Laws |
|---|---|
| Kernel | Lag range ↔ origin + tap count |
| Convolution | Delta identity; linearity; full commutativity; associativity (scale-aware tol); length/axis |
| Plans | Direct ≡ FFT ≡ overlap-add |
| Operator | Adjoint identity |
| Correlation | ≡ convolution with lag reversal |
| FIR runner | Batch from rest ≡ convolution |
| Filter runners | Chunking/concatenation; snapshot/restart |
| SOS | Cascade response ≡ product of sections |
| Representations | ZPK/SOS/TF/SS response agreement after conversion |
| FFT | Round-trip; Parseval; conjugate symmetry (real); circular conv theorem |
| STFT | Analysis/synthesis round-trip; axes match plan |
| Resampling | Ratio reduction; identity ratio ≡ identity (± delay) |
| Ownership | Immutable results stable after workspace reuse |
| Runners | Independent state per runner from one plan |
| Backends | Portable ≡ accelerated within contract |

Notes:

- Approximate floating-point laws with scale-sensitive tolerances.
- Do **not** publish Cats `Semiring` instances claiming exact Double associativity.
- Keep scalar reference implementations beside optimized kernels.
- Versioned SciPy/`gsignal` fixtures at generate-time; no Python/R at test runtime.
- External parity is evidence, not the API spec.

---

## 10. Implementation sequence

| PR | Focus |
|---|---|
| 1 | Build, core + laws, errors, rate/frequency, Sampling/Signal/Kernel, doc skeleton |
| 2 | Direct convolution (Full/Valid/Input), boundaries, axes, correlation, small-vector tests |
| 3 | Zero-boundary Gale operator + adjoint tests + direct benchmarks (no FFT) |
| 4 | FIR runner, Biquad, SOS (DF-IIt), snapshots, chunk laws, high-order fixtures |
| 5 | ComplexVector, FFT plans, mixed-radix, Bluestein, workspaces, JVM+JS benches |
| 6 | FFT convolution, OLA, measured Auto, method parity |
| 7 | Windows, FramePlan, StftPlan + dual, periodogram, Welch, typed spectra |
| 8 | FIR design, Butterworth, bilinear, ZPK→SOS, reports, response laws |
| 9 | RateRatio, polyphase, upfirdn, rational resampler, latency/chunk laws |
| 10 | Ravel adapters, API/numerical audit, benchmark receipts, migration examples, support matrix |

---

## 11. Review notes

### Strengths

- Clear separation of meaning vs method vs state matches how DSP actually fails in APIs (silent alignment, string modes, hidden filter state).
- Explicit kernel origin and `OutputRegion.Input` remove a real SciPy/`mode="same"` ambiguity.
- Gale operator + adjoint gives a path to inverse problems without inventing a second linear-algebra stack.
- Laws module and reference kernels make breadth growth safe.
- Deferred list is disciplined; 1.0 scope is large but architecturally coherent.

### Risks / watch points

1. **1.0 scope is still heavy** (FFT + STFT + design + multirate). Protect the first-slice sequence; resist starting at Butterworth or FFT.
2. **Boundary adjoints** for Reflect/Clamp are easy to ship half-correct—gate them behind tests.
3. **`Auto` convolution** needs real benchmark infrastructure early, or it becomes folklore anyway.
4. **STFT invertibility validation** is subtle; keep frame/window/hop validation in the plan constructor.
5. **Ravel optional** is right; do not let N-D batch convenience leak into core types.
6. **Analog/digital phantom types** are useful if kept tiny; avoid type-level sprawl beyond Domain tagging.

### Recommendation

Adopt the constitution and spine as written. Start PR 1–3 immediately. Treat `gsignal`/SciPy as fixture generators and coverage checklists only. Success metric for the early library: convolution regions, kernel origins, and the Gale adjoint law are boringly correct—then grow FFT, design, and multirate on that substrate.
