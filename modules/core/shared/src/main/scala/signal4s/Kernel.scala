package signal4s

import gale.linalg.DVec

/** Finite impulse-response kernel with an explicit zero-lag tap index.
  *
  * Tap `j` has lag `j - zeroLagIndex`. Even-length kernels must choose an
  * origin explicitly; [[Kernel.centeredOdd]] rejects even lengths.
  */
final case class Kernel private (
    taps: DVec,
    zeroLagIndex: Int
):
  def length: Int = taps.length
  def minLag: Int = -zeroLagIndex
  def maxLag: Int = taps.length - 1 - zeroLagIndex

  /** Time of the first full-convolution output sample relative to the input start. */
  def outputStartOffset(sampleRate: SampleRate): Seconds =
    Seconds.unsafe(-zeroLagIndex.toDouble / sampleRate.hertz)

object Kernel:
  def causal(taps: DVec): Either[SignalError, Kernel] =
    at(taps, zeroLagIndex = 0)

  def at(taps: DVec, zeroLagIndex: Int): Either[SignalError, Kernel] =
    if taps.length <= 0 then Left(SignalError.EmptyKernel)
    else if zeroLagIndex < 0 || zeroLagIndex >= taps.length then
      Left(SignalError.InvalidKernelOrigin(zeroLagIndex, taps.length))
    else Right(new Kernel(taps, zeroLagIndex))

  /** Center origin on the middle tap of an odd-length kernel. */
  def centeredOdd(taps: DVec): Either[SignalError, Kernel] =
    if taps.length <= 0 then Left(SignalError.EmptyKernel)
    else if taps.length % 2 == 0 then Left(SignalError.KernelNotOddLength(taps.length))
    else at(taps, zeroLagIndex = taps.length / 2)
