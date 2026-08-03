# signal4s

A compact calculus of sampled signals and systems for Scala 3 (JVM and Scala.js).

`signal4s` is not a translation of SciPy or gsignal. It organizes convolution,
filtering, Fourier analysis, spectral estimation, and resampling around a small
set of precise values:

**Signal + Kernel + System Representation + Plan + Runner + Typed Result Axes**

Start with the [documentation index](site-docs/README.md). The architecture,
scope, and implementation history remain in [proposal.md](proposal.md).

## Status

This checkout is an implementation snapshot, not a published release. The
current 1.0 implementation scope covers axes, convolution (direct/FFT/OLA),
STFT/Welch, FIR/Butterworth design, multirate (`upfirdn` / `ResamplePoly`),
optional `signal4s-ravel` axis adapters, Gale operators, and filter runners.
The [support reference](site-docs/reference/support.md) records the important
boundaries and deferred areas.

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

The build uses Gale through an explicit override, a sibling `../gale` checkout,
or a pinned git revision. Override with:

```bash
sbt -Dsignal4s.gale.build=/path/to/gale test
```

```bash
sbt coreJVM/test fftJVM/test designJVM/test ravelJVM/test lawsJVM/test
sbt coreJS/test fftJS/test designJS/test ravelJS/test lawsJS/test
sbt docs/tlSite
```

`docs/tlSite` compiles the public Markdown guides with mdoc and renders the
local Laika site. It does not publish anything. There are no verified Maven
coordinates for this snapshot yet; use the source build when working from a
checkout.

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

val sampleRate = SampleRate.hertz(1000.0).orThrow
val input = Signal(Vec(1.0, 2.0, 3.0, 4.0), Sampling(sampleRate)).orThrow
val kernel = Kernel.centeredOdd(Vec(0.25, 0.5, 0.25)).orThrow
val output = Convolution(input, kernel, OutputRegion.Full).orThrow

assert(output.length == 6)
assert(output.start == Seconds.of(-0.001).orThrow)
```

`Full` returns the complete finite-support convolution and therefore retains
the kernel's negative-lag coordinate. See the [getting-started guide](site-docs/getting-started.md)
for the region and boundary choices.
