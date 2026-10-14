package signal4s.fft

import gale.linalg.DVec
import gale.numeric.ExactSum
import signal4s.SignalError

/** These estimands differ by the weighted arithmetic/geometric baseline gap. */
enum LogPowerConvention:
  case LegacyLogRatio, MeanLogRatio

/** Power conventions only. Decibels use 10 log10, never amplitude's 20 log10. */
enum LogPowerScale:
  case NaturalLog, Log10, Decibels
  def base: Double = if this == NaturalLog then math.E else 10.0
  def multiplier: Double = if this == Decibels then 10.0 else 1.0
  private[fft] def log(value: Double): Double =
    if this == NaturalLog then math.log(value) else math.log10(value)

/** Owned scalar baseline; no time, frequency, trial or participant pooling is implicit. */
final class LogPowerBaseline private (val convention: LogPowerConvention,
    val scale: LogPowerScale, val baselineLogPower: Double,
    val selectedCount: Int, val positiveWeightCount: Int, val weightTotal: Double):
  def normalize(power: Double): Either[SignalError, Double] =
    if !power.isFinite || power <= 0 then Left(LogPowerBaseline.failure("strictly positive finite target power required"))
    else
      val result = scale.multiplier * (scale.log(power) - baselineLogPower)
      if result.isFinite then Right(result) else Left(LogPowerBaseline.failure("normalized log power exceeds finite capacity"))

object LogPowerBaseline:
  private[fft] def failure(detail: String): SignalError = SignalError.NumericalFailure("log-power-baseline", detail)

  /** Exact sums of once-rounded weighted terms. Masked powers are absent; zero
    * weights are selected but do not contribute. Admitted powers, including
    * zero-weight selected entries, must be strictly positive and finite.
    * Unrepresentable products/totals refuse; there is no flooring or rescaling.
    */
  def fit(powers: DVec, weights: DVec, valid: Vector[Boolean],
      convention: LogPowerConvention, scale: LogPowerScale): Either[SignalError, LogPowerBaseline] =
    if powers.length != weights.length || powers.length != valid.size then
      return Left(failure("power/weight/mask lengths differ"))
    val total = ExactSum.zero(); val normalizer = ExactSum.zero()
    var selected = 0; var positive = 0; var i = 0
    while i < powers.length do
      val weight = weights(i)
      if !weight.isFinite || weight < 0 then return Left(failure("finite nonnegative weights required"))
      if valid(i) then
        val power = powers(i)
        if !power.isFinite || power <= 0 then return Left(failure("strictly positive finite selected baseline power required"))
        selected += 1
        if weight > 0 then
          val value = convention match
            case LogPowerConvention.LegacyLogRatio => power
            case LogPowerConvention.MeanLogRatio => scale.log(power)
          val product = weight * value
          if !product.isFinite || (value != 0 && product == 0) then
            return Left(failure("weighted term overflows or underflows finite arithmetic capacity"))
          total.add(product) match
            case Left(error) => return Left(failure(error.message))
            case Right(_) => ()
          normalizer.add(weight) match
            case Left(error) => return Left(failure(error.message))
            case Right(_) => ()
          positive += 1
      i += 1
    if positive == 0 then Left(failure("at least one positive-weight admitted baseline observation required"))
    else if !normalizer.value.isFinite then Left(failure("reported weight total exceeds finite capacity"))
    else total.ratio(normalizer).left.map(e => failure(e.message)).flatMap: average =>
      val baseline = convention match
        case LogPowerConvention.LegacyLogRatio => scale.log(average)
        case LogPowerConvention.MeanLogRatio => average
      if !baseline.isFinite then Left(failure("baseline log power exceeds finite capacity"))
      else Right(new LogPowerBaseline(convention, scale, baseline, selected, positive, normalizer.value))
