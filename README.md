# signal4s

A compact calculus of sampled signals and systems for Scala 3 (JVM and Scala.js).

`signal4s` is not a translation of SciPy or gsignal. It organizes convolution,
filtering, Fourier analysis, spectral estimation, and resampling around a small
set of precise values:

**Signal + Kernel + System Representation + Plan + Runner + Typed Result Axes**

See [proposal.md](proposal.md) for the architecture, scope, and implementation
sequence.

## Status

**1.0 scope is implemented:** axes, convolution (direct/FFT/OLA), STFT/Welch,
FIR/Butterworth design, multirate (`upfirdn` / `ResamplePoly`), optional
`signal4s-ravel` axis adapters, Gale operator, and filter runners. SciPy/NumPy
parity uses committed fixtures under `fixtures/data/`. See
[docs/SUPPORT.md](docs/SUPPORT.md) and [docs/MIGRATION.md](docs/MIGRATION.md).

## Modules

| Module | Artifact | Role |
|---|---|---|
| `modules/core` | `signal4s-core` | Signals, axes, kernels, filters |
| `modules/fft` | `signal4s-fft` | FFT, convolution, windows, STFT, Welch |
| `modules/design` | `signal4s-design` | FIR/IIR design, bilinear, ZPK/SOS, responses |
| `modules/ravel` | `signal4s-ravel` | Optional multi-channel / along-axis adapters |
| `modules/laws` | `signal4s-laws` | Reusable munit/ScalaCheck law bundles |

Optional later: backends (Vector / native FFT), FS2 adapters.

## Build

Requires a Gale checkout. By default the build uses `../gale` when present,
or a pinned git revision. Override with:

```bash
sbt -Dsignal4s.gale.build=/path/to/gale test
```

```bash
sbt coreJVM/test fftJVM/test designJVM/test ravelJVM/test lawsJVM/test
sbt coreJS/test fftJS/test designJS/test ravelJS/test lawsJS/test
```

### Coverage

The JVM core and FFT quality gate requires at least 85% statement coverage.
Run it locally with:

```bash
sbt coverageCoreFft
```

The command runs the core and FFT JVM suites plus the reusable law suites, then
aggregates coverage for `signal4s-core` and `signal4s-fft`. Reports are written
to `target/scala-3.7.4/scoverage-report/`. Platform-specific JVM tests run on
macOS in CI, while the Scala.js core and FFT suites run separately.

## Performance vs SciPy

See [`docs/PERFORMANCE.md`](docs/PERFORMANCE.md). Comparative receipts live in
`benchmarks/receipts/vs-scipy.md`.

```bash
python3 fixtures/generate/bench_vs_scipy.py
sbt "fftJVM / Test / runMain signal4s.fft.bench.VsScipyBench"
```

## SciPy fixtures

Numeric parity against SciPy uses **versioned, committed fixtures**. Tests never
invoke Python. Regenerate offline from `fixtures/generate/` — see
[fixtures/README.md](fixtures/README.md). The JVM fixture loader validates every
committed smoke fixture, including its SciPy/NumPy pins and declared tolerances,
before a parity suite reads numeric values.

## Example

```scala
import gale.linalg.Vec
import signal4s.*

val fs = SampleRate.hertz(1000.0).orThrow
val x = Signal(
  samples = Vec.tabulate(1000)(i => math.sin(2.0 * math.Pi * 12.0 * i / fs.hertz)),
  sampling = Sampling(fs)
).orThrow

val smoother = Kernel.causal(Vec(0.25, 0.5, 0.25)).orThrow
val y = Convolution(x, smoother, OutputRegion.Full).orThrow
val _ = y
```
