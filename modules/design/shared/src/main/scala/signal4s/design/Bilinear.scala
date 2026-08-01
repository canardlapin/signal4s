package signal4s.design

import signal4s.{SampleRate, SignalError}
import signal4s.fft.Complex

/** Analog ↔ digital conversions via Tustin's method (SciPy `bilinear_zpk`). */
object Bilinear:

  /** Digital ZPK from analog ZPK at sample rate `fs` (no prewarping here).
    *
    * Substitutes \( s = 2 f_s (z-1)/(z+1) \). Relative degree zeros move to \( z=-1 \).
    */
  def transform(
      analog: ZerosPolesGain,
      sampleRate: SampleRate
  ): Either[SignalError, ZerosPolesGain] =
    val fs2 = 2.0 * sampleRate.hertz
    val zIn = analog.zeros
    val pIn = analog.poles
    val degree = pIn.length - zIn.length
    if degree < 0 then
      Left(
        SignalError.NumericalFailure(
          "Bilinear.transform",
          s"improper transfer function (zeros ${zIn.length} > poles ${pIn.length})"
        )
      )
    else
      val zOut = IArray.newBuilder[Complex]
      var i = 0
      while i < zIn.length do
        zOut += bilinearMap(zIn(i), fs2)
        i += 1
      var d = 0
      while d < degree do
        zOut += Complex(-1.0, 0.0)
        d += 1
      val pOut = IArray.newBuilder[Complex]
      i = 0
      while i < pIn.length do
        pOut += bilinearMap(pIn(i), fs2)
        i += 1
      val kNum = product(fs2, zIn)
      val kDen = product(fs2, pIn)
      if kDen.abs == 0.0 then
        Left(SignalError.NumericalFailure("Bilinear.transform", "pole at s=2*fs"))
      else
        val k = analog.gain * div(kNum, kDen).real
        ZerosPolesGain(zOut.result(), pOut.result(), k)

  private def bilinearMap(s: Complex, fs2: Double): Complex =
    // (fs2 + s) / (fs2 - s)
    val num = Complex(fs2 + s.real, s.imaginary)
    val den = Complex(fs2 - s.real, -s.imaginary)
    div(num, den)

  private def product(fs2: Double, roots: IArray[Complex]): Complex =
    var acc = Complex.One
    var i = 0
    while i < roots.length do
      acc = acc * Complex(fs2 - roots(i).real, -roots(i).imaginary)
      i += 1
    acc

  private def div(a: Complex, b: Complex): Complex =
    val n = b.real * b.real + b.imaginary * b.imaginary
    Complex(
      (a.real * b.real + a.imaginary * b.imaginary) / n,
      (a.imaginary * b.real - a.real * b.imaginary) / n
    )
