# Design a filter and estimate a spectrum

This is the common analysis path for a sampled measurement: design a stable
Butterworth filter as second-order sections, process the samples, then estimate
a Welch power spectrum.

```scala mdoc:silent
import gale.linalg.Vec
import signal4s.*
import signal4s.design.*
import signal4s.fft.*

val filterGuideRate = SampleRate.hertz(1000.0).orThrow
val filterGuideInput =
  Vec.tabulate(256)(i => math.sin(2 * math.Pi * 40 * i / filterGuideRate.hertz))

val filterGuideIir =
  Butterworth.lowPass(4, Frequency.hertz(100.0).orThrow, filterGuideRate).orThrow
val filterGuideOutput = filterGuideIir.sos.process(filterGuideInput).orThrow

val filterGuideWindow =
  Window.fromSpec(WindowSpec.Hann(64, WindowConvention.Periodic)).orThrow
val filterGuideWelch =
  WelchPlan(
    filterGuideWindow,
    hop = 32,
    sampleRate = filterGuideRate,
    nfft = 64,
    detrend = Detrend.Mean
  ).orThrow
val filterGuideSpectrum = filterGuideWelch.estimate(filterGuideOutput).orThrow

assert(filterGuideSpectrum.frequencies.length == 33)
assert(filterGuideSpectrum.power.length == 33)
```

The SOS representation avoids the numerical sensitivity of expanding a high-
order IIR into one polynomial transfer function. The window convention is
explicit: `Periodic` is the usual FFT-bin convention, while `Symmetric` is a
different finite-window definition.

The current Welch implementation is onesided. A `Twosided` plan can be
constructed but estimation returns a `NumericalFailure`; choose the default
`SpectralSides.Onesided` until that capability is implemented. The same
limitation applies to `Periodogram`.

For zero-phase filtering, use an implemented edge treatment such as
`EdgeTreatment.OddPad`. The enum also names SciPy's `Gustafsson` method, but
the current implementation does not provide a true Gustafsson calculation for
all filter forms; see [support and maturity](../reference/support.md).

Next: [carry a filter across chunks](streaming.md) or [reuse FFT plans](fft-and-plans.md).
