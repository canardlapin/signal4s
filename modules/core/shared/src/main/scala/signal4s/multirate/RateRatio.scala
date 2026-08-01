package signal4s.multirate

import signal4s.SignalError

/** Rational rate change \(f_{\mathrm{out}} / f_{\mathrm{in}} = \texttt{up}/\texttt{down}\).
  *
  * Always stored GCD-reduced with positive integers.
  */
final case class RateRatio private (up: Int, down: Int):
  def isIdentity: Boolean = up == 1 && down == 1
  def invert: RateRatio = new RateRatio(down, up)
  def scale: Double = up.toDouble / down.toDouble

object RateRatio:
  def apply(up: Int, down: Int): Either[SignalError, RateRatio] =
    if up <= 0 then Left(SignalError.InvalidInputLength(up))
    else if down <= 0 then Left(SignalError.InvalidInputLength(down))
    else
      val g = gcd(up, down)
      Right(new RateRatio(up / g, down / g))

  def identity: RateRatio = new RateRatio(1, 1)

  private[multirate] def unsafe(up: Int, down: Int): RateRatio =
    new RateRatio(up, down)

  private def gcd(a: Int, b: Int): Int =
    var x = a
    var y = b
    while y != 0 do
      val t = y
      y = x % y
      x = t
    x
