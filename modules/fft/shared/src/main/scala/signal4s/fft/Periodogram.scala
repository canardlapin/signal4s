package signal4s.fft

import gale.linalg.DVec
import signal4s.*

/** Single-segment periodogram (Welch with one frame covering the signal, or
  * windowed FFT of the full observation padded/truncated to `nfft`).
  */
object Periodogram:

  def apply(
      signal: DVec,
      sampleRate: SampleRate,
      window: Window,
      nfft: Int,
      detrend: Detrend = Detrend.Mean,
      sides: SpectralSides = SpectralSides.Onesided,
      scaling: SpectralScaling = SpectralScaling.Density
  ): Either[SignalError, WelchResult] =
    if signal.length != window.length then
      Left(
        SignalError.NumericalFailure(
          "Periodogram",
          s"window length ${window.length} must equal signal length ${signal.length} (SciPy periodogram pads window to nfft separately; signal is windowed at its length then zero-padded to nfft)"
        )
      )
    else
      // SciPy periodogram: window of len(x), FFT length nfft
      WelchPlan(
        window,
        hop = math.max(1, signal.length),
        sampleRate,
        nfft,
        detrend,
        sides,
        scaling,
        AverageMethod.Mean
      ).flatMap(_.estimate(signal))

  def hann(
      signal: DVec,
      sampleRate: SampleRate,
      nfft: Option[Int] = None,
      detrend: Detrend = Detrend.Mean,
      scaling: SpectralScaling = SpectralScaling.Density
  ): Either[SignalError, WelchResult] =
    val n = nfft.getOrElse(signal.length)
    Window
      .fromSpec(WindowSpec.Hann(signal.length, WindowConvention.Periodic))
      .flatMap(w => apply(signal, sampleRate, w, n, detrend, SpectralSides.Onesided, scaling))
