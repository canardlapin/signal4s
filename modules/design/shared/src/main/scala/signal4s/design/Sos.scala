package signal4s.design

import signal4s.SignalError
import signal4s.filter.{Biquad, SecondOrderCascade}
import signal4s.fft.Complex

/** Conversions into second-order sections. */
object Sos:

  /** Pair ZPK into biquads (SciPy-like conjugate / nearest pairing).
    *
    * Section order may differ from SciPy `zpk2sos`; frequency response should match.
    * Overall gain is placed in the first section.
    */
  def fromZerosPolesGain(zpk: ZerosPolesGain): Either[SignalError, SecondOrderCascade] =
    val poles = Array.tabulate(zpk.poles.length)(i => zpk.poles(i))
    val zeros = Array.tabulate(zpk.zeros.length)(i => zpk.zeros(i))
    val polePairs = pairRoots(poles)
    val zeroPairs = pairRoots(zeros)
    val nSec = math.max(polePairs.length, zeroPairs.length)
    if nSec == 0 then
      Biquad(zpk.gain, 0.0, 0.0, 0.0, 0.0).flatMap(bq => SecondOrderCascade.of(bq))
    else
      val sections = IArray.newBuilder[Biquad]
      var s = 0
      var err: Option[SignalError] = None
      while s < nSec && err.isEmpty do
        val (p1, p2) =
          if s < polePairs.length then polePairs(s)
          else (Complex.Zero, Complex.Zero)
        val (z1, z2) =
          if s < zeroPairs.length then zeroPairs(s)
          else (Complex.Zero, Complex.Zero)
        // Section gain: put overall gain on first section only
        val g = if s == 0 then zpk.gain else 1.0
        sectionFromPair(z1, z2, p1, p2, g) match
          case Left(e)  => err = Some(e)
          case Right(bq) => sections += bq
        s += 1
      err match
        case Some(e) => Left(e)
        case None    => SecondOrderCascade(1.0, sections.result())

  private def sectionFromPair(
      z1: Complex,
      z2: Complex,
      p1: Complex,
      p2: Complex,
      gain: Double
  ): Either[SignalError, Biquad] =
    // (z - z1)(z - z2) = z^2 - (z1+z2)z + z1 z2  → b coeffs with leading gain
    // poles → a coeffs; a0 = 1
    val b1raw = -(z1.real + z2.real) // imag parts cancel for conjugates / reals
    val b2raw = (z1 * z2).real
    val a1 = -(p1.real + p2.real)
    val a2 = (p1 * p2).real
    Biquad(gain, gain * b1raw, gain * b2raw, a1, a2)

  /** Pair conjugate roots; leftover reals paired with each other (pad with 0). */
  private def pairRoots(roots: Array[Complex]): Vector[(Complex, Complex)] =
    val used = Array.fill(roots.length)(false)
    val out = Vector.newBuilder[(Complex, Complex)]
    // First: conjugate pairs (positive imag first)
    var i = 0
    while i < roots.length do
      if !used(i) && math.abs(roots(i).imaginary) > 1e-12 then
        var j = i + 1
        var found = -1
        while j < roots.length && found < 0 do
          if !used(j) &&
            math.abs(roots(i).real - roots(j).real) < 1e-9 &&
            math.abs(roots(i).imaginary + roots(j).imaginary) < 1e-9
          then found = j
          j += 1
        if found >= 0 then
          val a = roots(i)
          val b = roots(found)
          if a.imaginary >= b.imaginary then out += ((a, b)) else out += ((b, a))
          used(i) = true
          used(found) = true
      i += 1
    // Remaining reals
    val reals = Vector.newBuilder[Complex]
    i = 0
    while i < roots.length do
      if !used(i) then reals += Complex(roots(i).real, 0.0)
      i += 1
    val rs = reals.result()
    var r = 0
    while r + 1 < rs.length do
      out += ((rs(r), rs(r + 1)))
      r += 2
    if r < rs.length then out += ((rs(r), Complex.Zero))
    out.result()
