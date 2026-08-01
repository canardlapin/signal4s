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
  private val MinDenom = 1e-12

  def fromAnalysis(analysis: Window, hop: Int): Either[SignalError, CanonicalDual] =
    if hop <= 0 then Left(SignalError.InvalidInputLength(hop))
    else
      val w = analysis.taps
      val m = w.length
      val denom = DVecBuilder.zeros(m)
      var n = 0
      while n < m do
        var s = 0.0
        var ell = -((m / hop) + 2)
        val ellMax = (m / hop) + 2
        while ell <= ellMax do
          val i = n + ell * hop
          if i >= 0 && i < m then s += w(i) * w(i)
          ell += 1
        denom(n) = s
        n += 1
      var ok = true
      var i = 0
      while i < m && ok do
        if denom(i) < MinDenom && math.abs(w(i)) > MinDenom then ok = false
        i += 1
      if !ok then
        Left(
          SignalError.NumericalFailure(
            "CanonicalDual",
            s"analysis window with hop=$hop is not invertible (vanishing overlap energy)"
          )
        )
      else
        val dual = DVecBuilder.zeros(m)
        i = 0
        while i < m do
          dual(i) = if denom(i) < MinDenom then 0.0 else w(i) / denom(i)
          i += 1
        Right(new CanonicalDual(dual.result(), hop))
