package signal4s.fft

import gale.linalg.DVec
import signal4s.SignalError
import signal4s.fft.internal.{FftEngine, Radix2Tables}

/** Immutable array layout. Inspection allocates no FFT buffers or tables. */
private[fft] final class OwnedFrameFft private (
    val nfft: Int, val workingLength: Int, val radix2: Boolean
):
  val workspaceBytes: BigInt =
    val n = BigInt(workingLength)
    val tables = 4 * n + 32 * (n - 1)
    val arrays = 16 * n + (if radix2 then BigInt(0) else 16 * n + 16 * BigInt(nfft))
    tables + arrays

  def newWorkspace(): OwnedFrameFftWorkspace = new OwnedFrameFftWorkspace(this)

object OwnedFrameFft:
  val MaxFftLength: Int = 1 << 28
  def apply(nfft: Int): Either[SignalError, OwnedFrameFft] =
    if nfft <= 0 || nfft > MaxFftLength then
      Left(SignalError.NumericalFailure("spectral.fft.shape", "FFT length outside bounded capacity"))
    else
      val radix2 = (nfft & (nfft - 1)) == 0
      var working = 1
      val minimum = if radix2 then nfft.toLong else 2L * nfft - 1
      while working.toLong < minimum do working <<= 1
      Right(new OwnedFrameFft(nfft, working, radix2))

/** Single-owner portable radix2/Bluestein state, shared by STFT and periodograms. */
private[fft] final class OwnedFrameFftWorkspace(val layout: OwnedFrameFft):
  private[fft] val re = Array.ofDim[Double](layout.workingLength)
  private[fft] val im = Array.ofDim[Double](layout.workingLength)
  private[fft] val tables = Radix2Tables.owned(layout.workingLength)
  private[fft] val chirpRe = if layout.radix2 then Array.emptyDoubleArray else Array.ofDim[Double](layout.nfft)
  private[fft] val chirpIm = if layout.radix2 then Array.emptyDoubleArray else Array.ofDim[Double](layout.nfft)
  private[fft] val kernelRe = if layout.radix2 then Array.emptyDoubleArray else Array.ofDim[Double](re.length)
  private[fft] val kernelIm = if layout.radix2 then Array.emptyDoubleArray else Array.ofDim[Double](re.length)
  if !layout.radix2 then
    var i = 0
    while i < layout.nfft do
      val angle = -math.Pi * ((i.toLong * i) % (2L * layout.nfft)).toDouble / layout.nfft.toDouble
      chirpRe(i) = math.cos(angle)
      chirpIm(i) = math.sin(angle)
      kernelRe(i) = chirpRe(i)
      kernelIm(i) = -chirpIm(i)
      if i > 0 then
        kernelRe(re.length - i) = chirpRe(i)
        kernelIm(re.length - i) = -chirpIm(i)
      i += 1
    FftEngine.radix2Owned(kernelRe, kernelIm, false, tables)

  def forward(frame: DVec): Unit =
    require(frame.length <= layout.nfft, "prepared frame exceeds FFT length")
    java.util.Arrays.fill(re, 0.0)
    java.util.Arrays.fill(im, 0.0)
    var i = 0
    while i < frame.length do
      if layout.radix2 then re(i) = frame(i)
      else
        re(i) = frame(i) * chirpRe(i)
        im(i) = frame(i) * chirpIm(i)
      i += 1
    FftEngine.radix2Owned(re, im, false, tables)
    if !layout.radix2 then
      i = 0
      while i < re.length do
        val ar = re(i)
        val ai = im(i)
        re(i) = ar * kernelRe(i) - ai * kernelIm(i)
        im(i) = ar * kernelIm(i) + ai * kernelRe(i)
        i += 1
      FftEngine.radix2Owned(re, im, true, tables)

  def real(k: Int): Double =
    if layout.radix2 then re(k)
    else (re(k) * chirpRe(k) - im(k) * chirpIm(k)) / re.length

  def imaginary(k: Int): Double =
    if layout.radix2 then im(k)
    else (re(k) * chirpIm(k) + im(k) * chirpRe(k)) / re.length
