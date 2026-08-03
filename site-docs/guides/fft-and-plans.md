# FFT backends and reusable plans

The FFT module installs optional convolution hooks into `signal4s-core`. Make
that installation explicit at application startup when using FFT or overlap-add
convolution through the core API.

```scala mdoc:silent
import gale.linalg.Vec
import signal4s.*
import signal4s.fft.*

FftBackend.ensureInstalled()

val fftGuideRate = SampleRate.hertz(1000.0).orThrow
val fftGuideRealPlan =
  RealFftPlan(8, FftNormalization.Backward, fftGuideRate).orThrow
val fftGuideSpectrum =
  fftGuideRealPlan.forward(Vec(0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0)).orThrow
assert(fftGuideSpectrum.bins.length == 5)

val fftGuideKernel = Kernel.causal(Vec(0.25, 0.5, 0.25)).orThrow
val fftGuideConvolutionPlan =
  Convolution
    .plan(
      fftGuideKernel,
      inputLength = 1024,
      region = OutputRegion.Full,
      method = ConvolutionMethod.Fft
    )
    .orThrow
val fftGuideOutput =
  fftGuideConvolutionPlan.apply(Vec.tabulate(1024)(_ => 1.0)).orThrow
assert(fftGuideOutput.length == 1026)
```

`FftPlan` and `RealFftPlan` are immutable and fixed to a transform length.
Their `newWorkspace()` values are mutable scratch and must remain single-owner.
The same distinction applies to FFT-backed `ConvolutionPlan` values: the plan
captures the kernel, input length, region, and selected method, while execution
receives the sample vector.

Importing `signal4s.fft.*` alone is not a reliable registration contract for
core `Convolution(..., ConvolutionMethod.Auto)`. Call
`FftBackend.ensureInstalled()` explicitly before constructing or invoking an
FFT-backed core operation. If hooks are absent, core can run `Direct` for
`Auto` or returns an error for an explicitly requested FFT method.

Next: [interpret the support boundaries](../reference/support.md) and the
[scoped performance evidence](../reference/performance.md).
