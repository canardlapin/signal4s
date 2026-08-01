package signal4s.design

import signal4s.SignalError
import signal4s.fft.Complex

/** Zeros, poles, and gain of a rational transfer function. */
final case class ZerosPolesGain private (
    zeros: IArray[Complex],
    poles: IArray[Complex],
    gain: Double
):
  def order: Int = math.max(zeros.length, poles.length)

object ZerosPolesGain:
  def apply(
      zeros: IArray[Complex],
      poles: IArray[Complex],
      gain: Double
  ): Either[SignalError, ZerosPolesGain] =
    if !gain.isFinite then Left(SignalError.InvalidFilterCoefficient("gain", gain))
    else if zeros.exists(z => !z.real.isFinite || !z.imaginary.isFinite) then
      Left(SignalError.NumericalFailure("ZerosPolesGain", "non-finite zero"))
    else if poles.exists(p => !p.real.isFinite || !p.imaginary.isFinite) then
      Left(SignalError.NumericalFailure("ZerosPolesGain", "non-finite pole"))
    else Right(new ZerosPolesGain(zeros, poles, gain))

  def unsafe(zeros: IArray[Complex], poles: IArray[Complex], gain: Double): ZerosPolesGain =
    new ZerosPolesGain(zeros, poles, gain)
