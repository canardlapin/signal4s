package signal4s.fft

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*
import signal4s.fft.internal.{FastFftLength, FftEngine}

/** Linear convolution via FFT (zero-padded circular convolution). */
private[fft] object FftConvolution:

  /** Next transform length ≥ `min` (power-of-two or 5-smooth; see [[FastFftLength]]). */
  def fftLength(min: Int): Int =
    FastFftLength(math.max(2, min))

  private def isPowerOfTwo(n: Int): Boolean = FastFftLength.isPowerOfTwo(n)

  def full(signal: DVec, kernel: Kernel): Either[SignalError, DVec] =
    val outLen = signal.length + kernel.length - 1
    val nfft = fftLength(outLen)
    convolveFullTo(signal, kernel, nfft).map(y => take(y, outLen))

  /** Apply a pre-transformed kernel spectrum (length `nfft`). */
  def fullWithKernelSpectrum(
      signal: DVec,
      kernelSpectrum: ComplexVector,
      outLen: Int
  ): Either[SignalError, DVec] =
    val nfft = kernelSpectrum.length
    val kr = new Array[Double](nfft)
    val ki = new Array[Double](nfft)
    kernelSpectrum.real.copyTo(kr)
    kernelSpectrum.imaginary.copyTo(ki)
    val (hr, hi) = onesidedFromHermitian(kr, ki)
    fullWithKernelSpectrumArrays(signal, hr, hi, outLen)

  /** Like [[fullWithKernelSpectrumArrays]] but reuses caller-owned FFT work buffers.
    *
    * `workRe` / `workIm` must each have length `kernelRe.length`. `packRe` /
    * `packIm` must each have length `nfft/2` for the packed real FFT. Not
    * thread-safe with respect to those buffers.
    */
  def fullWithKernelSpectrumArraysInto(
      signal: DVec,
      kernelHalfRe: Array[Double],
      kernelHalfIm: Array[Double],
      outLen: Int,
      workRe: Array[Double],
      workIm: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double]
  ): Either[SignalError, DVec] =
    fullWithKernelSpectrumArraysInto(
      signal,
      kernelHalfRe,
      kernelHalfIm,
      outLen,
      workRe,
      workIm,
      packRe,
      packIm,
      new Array[Double](outLen)
    )

  def fullWithKernelSpectrumArraysInto(
      signal: DVec,
      kernelHalfRe: Array[Double],
      kernelHalfIm: Array[Double],
      outLen: Int,
      workRe: Array[Double],
      workIm: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double],
      outBuf: Array[Double]
  ): Either[SignalError, DVec] =
    fullBlockWithKernelSpectrumArraysInto(
      signal,
      signalOffset = 0,
      signalLength = signal.length,
      kernelHalfRe,
      kernelHalfIm,
      outLen,
      workRe,
      workIm,
      packRe,
      packIm,
      outBuf
    ).map(_ => adopt(outBuf, outLen))

  /** Transform one contiguous signal block into caller-owned output and FFT buffers. */
  private[fft] def fullBlockWithKernelSpectrumArraysInto(
      signal: DVec,
      signalOffset: Int,
      signalLength: Int,
      kernelHalfRe: Array[Double],
      kernelHalfIm: Array[Double],
      outLen: Int,
      workRe: Array[Double],
      workIm: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double],
      outBuf: Array[Double]
  ): Either[SignalError, Unit] =
    val bins = kernelHalfRe.length
    if kernelHalfIm.length != bins || bins < 2 then
      Left(SignalError.LengthMismatch(bins, kernelHalfIm.length))
    else
      val nfft = (bins - 1) << 1
      if workRe.length != nfft || workIm.length != nfft then
        Left(SignalError.LengthMismatch(nfft, workRe.length))
      else if outBuf.length < outLen then
        Left(SignalError.LengthMismatch(outLen, outBuf.length))
      else if outLen < 0 || outLen > nfft then
        Left(SignalError.NumericalFailure("FftConvolution", s"outLen=$outLen nfft=$nfft"))
      else if signalOffset < 0 || signalLength < 0 || signalOffset + signalLength > signal.length then
        Left(
          SignalError.NumericalFailure(
            "FftConvolution",
            s"invalid block offset=$signalOffset length=$signalLength signalLength=${signal.length}"
          )
        )
      else if signalLength > nfft then
        Left(SignalError.LengthMismatch(nfft, signalLength))
      else if isPowerOfTwo(nfft) then
        if packRe.length != nfft / 2 || packIm.length != nfft / 2 then
          Left(SignalError.LengthMismatch(nfft / 2, packRe.length))
        else
          copyBlock(signal, signalOffset, signalLength, workRe)
          if signalLength < nfft then
            java.util.Arrays.fill(workRe, signalLength, nfft, 0.0)
          // Onesided spectrum: inverse packed real path only reads bins 0..nfft/2.
          FftEngine.forwardRealOnesided(workRe, workIm, packRe, packIm)
          multiplyOnesided(workRe, workIm, kernelHalfRe, kernelHalfIm)
          FftEngine.inverseRealScaledTo(
            workRe,
            workIm,
            packRe,
            packIm,
            scale = 1.0 / nfft.toDouble,
            outBuf,
            outLen
          )
          Right(())
      else
        // 5-smooth (or other engine-supported) length: complex FFT of a real signal.
        copyBlock(signal, signalOffset, signalLength, workRe)
        if signalLength < nfft then
          java.util.Arrays.fill(workRe, signalLength, nfft, 0.0)
        java.util.Arrays.fill(workIm, 0.0)
        FftEngine.forward(workRe, workIm)
        multiplyHermitianHalf(workRe, workIm, kernelHalfRe, kernelHalfIm)
        FftEngine.inverse(workRe, workIm)
        val scale = 1.0 / nfft.toDouble
        var i = 0
        while i < outLen do
          outBuf(i) = workRe(i) * scale
          i += 1
        Right(())

  def fullWithKernelSpectrumArrays(
      signal: DVec,
      kernelHalfRe: Array[Double],
      kernelHalfIm: Array[Double],
      outLen: Int
  ): Either[SignalError, DVec] =
    val nfft = (kernelHalfRe.length - 1) << 1
    val packN = if isPowerOfTwo(nfft) then nfft / 2 else 0
    fullWithKernelSpectrumArraysInto(
      signal,
      kernelHalfRe,
      kernelHalfIm,
      outLen,
      new Array[Double](nfft),
      new Array[Double](nfft),
      new Array[Double](packN),
      new Array[Double](packN)
    )

  def transformKernel(kernel: Kernel, nfft: Int): Either[SignalError, ComplexVector] =
    if nfft <= 0 then Left(SignalError.InvalidInputLength(nfft))
    else if kernel.length > nfft then Left(SignalError.LengthMismatch(nfft, kernel.length))
    else
      val re = new Array[Double](nfft)
      val im = new Array[Double](nfft)
      kernel.taps.copyTo(re)
      if isPowerOfTwo(nfft) then
        val packRe = new Array[Double](nfft / 2)
        val packIm = new Array[Double](nfft / 2)
        FftEngine.forwardReal(re, im, packRe, packIm)
      else
        java.util.Arrays.fill(im, 0.0)
        FftEngine.forward(re, im)
      val reB = DVecBuilder.zeros(nfft)
      val imB = DVecBuilder.zeros(nfft)
      var i = 0
      while i < nfft do
        reB(i) = re(i)
        imB(i) = im(i)
        i += 1
      Right(ComplexVector.unsafe(reB.result(), imB.result()))

  /** Transform a real kernel directly to its Nyquist-inclusive onesided spectrum. */
  def transformKernelHalfArrays(
      kernel: Kernel,
      nfft: Int
  ): Either[SignalError, (Array[Double], Array[Double])] =
    if nfft <= 0 then Left(SignalError.InvalidInputLength(nfft))
    else if kernel.length > nfft then Left(SignalError.LengthMismatch(nfft, kernel.length))
    else
      val re = new Array[Double](nfft)
      val im = new Array[Double](nfft)
      kernel.taps.copyTo(re)
      if isPowerOfTwo(nfft) then
        val packRe = new Array[Double](nfft / 2)
        val packIm = new Array[Double](nfft / 2)
        FftEngine.forwardRealOnesided(re, im, packRe, packIm)
      else FftEngine.forward(re, im)
      Right(onesidedFromHermitian(re, im))

  private def convolveFullTo(
      signal: DVec,
      kernel: Kernel,
      nfft: Int
  ): Either[SignalError, DVec] =
    for
      half <- transformKernelHalfArrays(kernel, nfft)
      y <- fullWithKernelSpectrumArrays(signal, half._1, half._2, signal.length + kernel.length - 1)
    yield y

  /** Pointwise product on onesided bins `0 .. n/2` (packed inverse reads only these). */
  private def multiplyOnesided(
      aRe: Array[Double],
      aIm: Array[Double],
      bHalfRe: Array[Double],
      bHalfIm: Array[Double]
  ): Unit =
    val n = aRe.length
    val N = n >> 1
    require(bHalfRe.length == N + 1 && bHalfIm.length == N + 1)
    aRe(0) *= bHalfRe(0)
    aIm(0) = 0.0
    aRe(N) *= bHalfRe(N)
    aIm(N) = 0.0
    var k = 1
    while k < N do
      val ar = aRe(k)
      val ai = aIm(k)
      val br = bHalfRe(k)
      val bi = bHalfIm(k)
      aRe(k) = ar * br - ai * bi
      aIm(k) = ar * bi + ai * br
      k += 1

  /** Full Hermitian product (complex / non-pow2 FFT convolution path). */
  private def multiplyHermitianHalf(
      aRe: Array[Double],
      aIm: Array[Double],
      bHalfRe: Array[Double],
      bHalfIm: Array[Double]
  ): Unit =
    val n = aRe.length
    val N = n >> 1
    require(bHalfRe.length == N + 1 && bHalfIm.length == N + 1)
    aRe(0) *= bHalfRe(0)
    aIm(0) = 0.0
    aRe(N) *= bHalfRe(N)
    aIm(N) = 0.0
    var k = 1
    while k < N do
      val ar = aRe(k)
      val ai = aIm(k)
      val br = bHalfRe(k)
      val bi = bHalfIm(k)
      val pr = ar * br - ai * bi
      val pi = ar * bi + ai * br
      aRe(k) = pr
      aIm(k) = pi
      aRe(n - k) = pr
      aIm(n - k) = -pi
      k += 1

  /** Onesided (rfft) view of a Hermitian spectrum: length `n/2+1`. */
  def onesidedFromHermitian(re: Array[Double], im: Array[Double]): (Array[Double], Array[Double]) =
    val N = re.length >> 1
    val hr = new Array[Double](N + 1)
    val hi = new Array[Double](N + 1)
    var k = 0
    while k <= N do
      hr(k) = re(k)
      hi(k) = im(k)
      k += 1
    (hr, hi)

  private def adopt(values: Array[Double], length: Int): DVec =
    val builder = DVecBuilder.zeros(length)
    var i = 0
    while i < length do
      builder(i) = values(i)
      i += 1
    builder.result()

  private def copyBlock(
      signal: DVec,
      offset: Int,
      length: Int,
      target: Array[Double]
  ): Unit =
    var i = 0
    while i < length do
      target(i) = signal(offset + i)
      i += 1

  private def take(x: DVec, n: Int): DVec =
    if x.length == n then x
    else
      val out = DVecBuilder.zeros(n)
      var i = 0
      while i < n do
        out(i) = x(i)
        i += 1
      out.result()

  /** Region extraction from a full convolution result (same meaning as direct). */
  def extractRegion(
      full: DVec,
      signalLength: Int,
      kernel: Kernel,
      region: OutputRegion
  ): Either[SignalError, DVec] =
    region match
      case OutputRegion.Full =>
        Right(full)
      case OutputRegion.Valid =>
        // Match [[signal4s.internal.DirectConvolution.valid]]: first sample index
        // is `m - 1 - z` on the origin-independent full array `full[i]=Σ h[j] x[i-j]`.
        val z = kernel.zeroLagIndex
        val first = kernel.length - 1 - z
        val last = signalLength - 1 - z
        val len = math.max(0, last - first + 1)
        if len == 0 then Right(DVec.zeros(0))
        else Right(full.slice(first, first + len).copy)
      case OutputRegion.Input(Boundary.Zero) =>
        // SciPy-style full[i]=Σ h[j]x[i-j]; Input_z[n]=Σ h[j]x[n-j+z]=full[n+z].
        val start = kernel.zeroLagIndex
        Right(full.slice(start, start + signalLength).copy)
      case OutputRegion.Input(other) =>
        Left(
          SignalError.NumericalFailure(
            "FftConvolution",
            s"FFT convolution supports Full, Valid, and Input(Zero); got Input($other)"
          )
        )
