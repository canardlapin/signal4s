package signal4s.laws

import gale.linalg.DVec
import signal4s.*
import signal4s.fft.FftBackend
import munit.Assertions

object MethodParityLaws extends Assertions:

  FftBackend.ensureInstalled()

  def directEqualsFft(
      signal: DVec,
      kernel: Kernel,
      region: OutputRegion,
      tol: Double = 1e-9
  ): Unit =
    val d = Convolution(signal, kernel, region, ConvolutionMethod.Direct).orThrow
    val f = Convolution(signal, kernel, region, ConvolutionMethod.Fft).orThrow
    ConvolutionLaws.assertClose(f, d, tol)

  def directEqualsOla(
      signal: DVec,
      kernel: Kernel,
      block: Int,
      region: OutputRegion,
      tol: Double = 1e-9
  ): Unit =
    val d = Convolution(signal, kernel, region, ConvolutionMethod.Direct).orThrow
    val o =
      Convolution(signal, kernel, region, ConvolutionMethod.OverlapAdd(block)).orThrow
    ConvolutionLaws.assertClose(o, d, tol)

  def planReusesKernel(
      signal: DVec,
      kernel: Kernel,
      region: OutputRegion = OutputRegion.Full
  ): Unit =
    val plan =
      Convolution.plan(kernel, signal.length, region, ConvolutionMethod.Fft).orThrow
    assertEquals(plan.selectedMethod, ConvolutionMethod.Fft)
    val y1 = plan(signal).orThrow
    val y2 = plan(signal).orThrow
    ConvolutionLaws.assertClose(y1, y2, 1e-12)
    val direct = Convolution(signal, kernel, region, ConvolutionMethod.Direct).orThrow
    ConvolutionLaws.assertClose(y1, direct, 1e-9)

  /** Immutable result remains stable after another plan application. */
  def ownershipStable(signal: DVec, kernel: Kernel): Unit =
    val plan =
      Convolution.plan(kernel, signal.length, OutputRegion.Full, ConvolutionMethod.Fft).orThrow
    val y1 = plan(signal).orThrow
    val snapshot = y1.copy
    val _ = plan(signal).orThrow
    ConvolutionLaws.assertClose(y1, snapshot, 0.0)

  def circularConvolutionTheorem(signal: DVec, kernel: Kernel): Unit =
    val n = signal.length
    val direct = Convolution.circular(signal, kernel, n).orThrow
    val viaFft = FftBackend.circular(signal, kernel, n).orThrow
    ConvolutionLaws.assertClose(viaFft, direct, 1e-9)
