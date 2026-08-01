package signal4s.fft

import gale.linalg.DVec
import signal4s.*
import signal4s.internal.{ConvolutionHooks, DirectConvolution}

/** Installs FFT/OLA convolution hooks into `signal4s-core`.
  *
  * Call [[ensureInstalled]] from any code path that needs FFT convolution
  * (laws and fft module tests do this automatically via class loading).
  */
object FftBackend:

  ensureInstalled()

  def ensureInstalled(): Unit =
    ConvolutionHooks.install(run, plan, runCircular)

  private def run(
      signal: DVec,
      kernel: Kernel,
      region: OutputRegion,
      method: ConvolutionMethod
  ): Either[SignalError, DVec] =
    PlannedConvolution.make(kernel, signal.length, region, method).flatMap(_.apply(signal))

  private def plan(
      kernel: Kernel,
      inputLength: Int,
      region: OutputRegion,
      method: ConvolutionMethod
  ): Either[SignalError, ConvolutionPlan] =
    PlannedConvolution.make(kernel, inputLength, region, method)

  /** Circular convolution via DFT product (for theorem laws / cross-checks). */
  def circular(signal: DVec, kernel: Kernel, period: Int): Either[SignalError, DVec] =
    CircularFft.convolve(signal, kernel, period)

  private def runCircular(signal: DVec, kernel: Kernel, period: Int): Either[SignalError, DVec] =
    if CircularFft.shouldUse(period, kernel.length) then CircularFft.convolve(signal, kernel, period)
    else Right(DirectConvolution.circular(signal, kernel, period))
