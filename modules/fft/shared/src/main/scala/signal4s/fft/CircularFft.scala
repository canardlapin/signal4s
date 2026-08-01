package signal4s.fft

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*
import signal4s.fft.internal.FftEngine

/** Circular convolution via the DFT convolution theorem. */
private[fft] object CircularFft:

  def convolve(signal: DVec, kernel: Kernel, period: Int): Either[SignalError, DVec] =
    if period <= 0 then Left(SignalError.InvalidPeriod(period))
    else if signal.length != period then Left(SignalError.LengthMismatch(period, signal.length))
    else if FftEngine.supportsRealLength(period) then
      convolveReal(signal, kernel, period)
    else
      for
        plan <- FftPlan(period, FftNormalization.Backward)
        xSpec <- forwardReal(signal, plan)
        hSpec <- forwardKernel(kernel, period, plan)
        prod <- mul(xSpec, hSpec)
        y <- plan.inverse(prod)
      yield y.real

  /** Conservative circular crossover for the installed backend.
    *
    * Native packed real transforms win across the measured 64–8192 period grid
    * once a nontrivial FIR is present. The portable path remains direct until
    * its \(O(NM)\) work exceeds an FFT-scale bound.
    */
  def shouldUse(period: Int, kernelLength: Int): Boolean =
    if !FftEngine.supportsRealLength(period) then false
    else if FftEngine.hasNativeRealBackend then period >= 64 && kernelLength >= 8
    else
      val logN = math.max(1, 32 - Integer.numberOfLeadingZeros(period - 1))
      kernelLength >= 8 * logN

  /** Real circular convolution avoids split-complex vectors and full spectra. */
  private def convolveReal(
      signal: DVec,
      kernel: Kernel,
      period: Int
  ): Either[SignalError, DVec] =
    val workRe = new Array[Double](period)
    val workIm = new Array[Double](period)
    val packRe = new Array[Double](period / 2)
    val packIm = new Array[Double](period / 2)
    signal.copyTo(workRe)
    FftEngine.forwardRealOnesided(workRe, workIm, packRe, packIm)

    val kernelRe = new Array[Double](period)
    val kernelIm = new Array[Double](period)
    val kernelPackRe = new Array[Double](period / 2)
    val kernelPackIm = new Array[Double](period / 2)
    val z = kernel.zeroLagIndex
    var j = 0
    while j < kernel.length do
      val idx = Math.floorMod(j - z, period)
      kernelRe(idx) += kernel.taps(j)
      j += 1
    FftEngine.forwardRealOnesided(kernelRe, kernelIm, kernelPackRe, kernelPackIm)

    multiplyOnesided(workRe, workIm, kernelRe, kernelIm)
    val out = new Array[Double](period)
    FftEngine.inverseRealScaledTo(
      workRe,
      workIm,
      packRe,
      packIm,
      scale = 1.0 / period.toDouble,
      out,
      period
    )
    val result = DVecBuilder.zeros(period)
    var i = 0
    while i < period do
      result(i) = out(i)
      i += 1
    Right(result.result())

  private def forwardReal(x: DVec, plan: FftPlan): Either[SignalError, ComplexVector] =
    val re = DVecBuilder.zeros(plan.length)
    val im = DVecBuilder.zeros(plan.length)
    var i = 0
    while i < x.length do
      re(i) = x(i)
      i += 1
    plan.forward(ComplexVector.unsafe(re.result(), im.result()))

  private def forwardKernel(
      kernel: Kernel,
      period: Int,
      plan: FftPlan
  ): Either[SignalError, ComplexVector] =
    // Place taps so circular convolution matches DirectConvolution.circular:
    // y[i] = Σ_j h[j] x[(i - (j-z)) mod N]
    val z = kernel.zeroLagIndex
    val re = DVecBuilder.zeros(period)
    val im = DVecBuilder.zeros(period)
    var j = 0
    while j < kernel.length do
      val lag = j - z
      val idx = Math.floorMod(lag, period)
      re(idx) = re(idx) + kernel.taps(j)
      j += 1
    plan.forward(ComplexVector.unsafe(re.result(), im.result()))

  private def mul(a: ComplexVector, b: ComplexVector): Either[SignalError, ComplexVector] =
    if a.length != b.length then Left(SignalError.LengthMismatch(a.length, b.length))
    else
      val re = DVecBuilder.zeros(a.length)
      val im = DVecBuilder.zeros(a.length)
      var i = 0
      while i < a.length do
        re(i) = a.real(i) * b.real(i) - a.imaginary(i) * b.imaginary(i)
        im(i) = a.real(i) * b.imaginary(i) + a.imaginary(i) * b.real(i)
        i += 1
      Right(ComplexVector.unsafe(re.result(), im.result()))

  private def multiplyOnesided(
      signalRe: Array[Double],
      signalIm: Array[Double],
      kernelRe: Array[Double],
      kernelIm: Array[Double]
  ): Unit =
    val half = signalRe.length / 2
    signalRe(0) *= kernelRe(0)
    signalIm(0) = 0.0
    signalRe(half) *= kernelRe(half)
    signalIm(half) = 0.0
    var k = 1
    while k < half do
      val ar = signalRe(k)
      val ai = signalIm(k)
      signalRe(k) = ar * kernelRe(k) - ai * kernelIm(k)
      signalIm(k) = ar * kernelIm(k) + ai * kernelRe(k)
      k += 1
