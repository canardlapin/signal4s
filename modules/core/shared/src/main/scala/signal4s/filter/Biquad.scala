package signal4s.filter

import signal4s.SignalError

/** Normalized second-order section (\(a_0 = 1\)).
  *
  * Difference equation:
  * \( y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] - a_1 y[n-1] - a_2 y[n-2] \).
  */
final case class Biquad private (
    b0: Double,
    b1: Double,
    b2: Double,
    a1: Double,
    a2: Double
):
  /** Direct-form II transposed state length (two delays). */
  def stateLength: Int = 2

object Biquad:
  def apply(
      b0: Double,
      b1: Double,
      b2: Double,
      a1: Double,
      a2: Double
  ): Either[SignalError, Biquad] =
    def finite(name: String, value: Double): Either[SignalError, Unit] =
      if value.isFinite then Right(())
      else Left(SignalError.InvalidFilterCoefficient(name, value))

    for
      _ <- finite("b0", b0)
      _ <- finite("b1", b1)
      _ <- finite("b2", b2)
      _ <- finite("a1", a1)
      _ <- finite("a2", a2)
    yield new Biquad(b0, b1, b2, a1, a2)

  /** Construct from SciPy-style row `[b0, b1, b2, a0, a1, a2]`, normalizing by `a0`. */
  def fromSosRow(row: IArray[Double]): Either[SignalError, Biquad] =
    if row.length != 6 then
      Left(SignalError.LengthMismatch(6, row.length))
    else
      val a0 = row(3)
      if !a0.isFinite || a0 == 0.0 then
        Left(SignalError.InvalidFilterCoefficient("a0", a0))
      else
        apply(
          b0 = row(0) / a0,
          b1 = row(1) / a0,
          b2 = row(2) / a0,
          a1 = row(4) / a0,
          a2 = row(5) / a0
        )
