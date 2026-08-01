package signal4s.design

import gale.linalg.DVec
import signal4s.SignalError
import signal4s.filter.DigitalTransferFunction
import signal4s.fft.Complex

/** Polynomial / TF conversions from ZPK. */
object Transfer:

  def fromZerosPolesGain(zpk: ZerosPolesGain): Either[SignalError, DigitalTransferFunction] =
    val bPoly = polyFromRoots(zpk.zeros).map(_ * zpk.gain)
    val aPoly = polyFromRoots(zpk.poles)
    // Highest power first → negative powers of z (b0 + b1 z^{-1} + …)
    val b = DVec.fromSeq(bPoly.reverse.toIndexedSeq)
    val a = DVec.fromSeq(aPoly.reverse.toIndexedSeq)
    DigitalTransferFunction(b, a)

  /** Monic polynomial coefficients, highest power first. */
  private def polyFromRoots(roots: IArray[Complex]): Array[Double] =
    var coeffs = Array(Complex.One)
    var r = 0
    while r < roots.length do
      val root = roots(r)
      val grown = Array.fill(coeffs.length + 1)(Complex.Zero)
      var i = 0
      while i < coeffs.length do
        grown(i) = grown(i) + coeffs(i)
        grown(i + 1) = grown(i + 1) - coeffs(i) * root
        i += 1
      coeffs = grown
      r += 1
    coeffs.map(_.real)
