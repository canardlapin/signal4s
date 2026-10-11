package signal4s.fft

import gale.linalg.{DVec, DVecBuilder}
import signal4s.SignalError

/** Canonical dual window for STFT overlap-add reconstruction.
  *
  * \(\gamma[n] = w[n] / \sum_\ell w[n + \ell H]^2\) for indices in range,
  * matching SciPy `ShortTimeFFT.dual_win`.
  */
final case class CanonicalDual private (taps: DVec, hop: Int):
  def length: Int = taps.length

object CanonicalDual:
  def fromAnalysis(analysis: Window, hop: Int): Either[SignalError, CanonicalDual] =
    if hop <= 0 then Left(SignalError.InvalidInputLength(hop))
    else
      val w = analysis.taps
      val m = w.length
      if hop > m then
        return Left(SignalError.NumericalFailure("CanonicalDual", "hop exceeds window length; uncovered overlap phases"))
      // Every hop residue must be covered, including zero taps. Scale each residue
      // separately so small/large valid windows do not depend on an absolute cutoff.
      val scales = Array.fill(hop)(0.0)
      var i = 0
      while i < m do
        val residue = i % hop
        scales(residue) = math.max(scales(residue), math.abs(w(i)))
        i += 1
      if scales.exists(_ == 0.0) then
        Left(SignalError.NumericalFailure("CanonicalDual", s"analysis window with hop=$hop has an uncovered overlap phase"))
      else
        val energy = Array.fill(hop)(0.0)
        i = 0
        while i < m do
          val residue = i % hop
          val normalized = w(i) / scales(residue)
          energy(residue) += normalized * normalized
          i += 1
        val dual = DVecBuilder.zeros(m)
        i = 0
        while i < m do
          val residue = i % hop
          val value = (w(i) / scales(residue)) / energy(residue) / scales(residue)
          if !value.isFinite then
            return Left(SignalError.NumericalFailure("CanonicalDual", "dual coefficient exceeds finite capacity"))
          dual(i) = value
          i += 1
        Right(new CanonicalDual(dual.result(), hop))
