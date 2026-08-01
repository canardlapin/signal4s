package signal4s

import gale.linalg.{DVec, DVecBuilder}

/** Output lag coordinates for a correlation. Lag at index `i` is `first + i`. */
final case class LagAxis private (first: Int, length: Int):
  def lagAt(index: Int): Int =
    require(index >= 0 && index < length, s"lag index $index out of range [0, $length)")
    first + index

  def last: Int = first + length - 1

object LagAxis:
  def apply(first: Int, length: Int): Either[SignalError, LagAxis] =
    if length < 0 then Left(SignalError.LengthMismatch(0, length))
    else Right(new LagAxis(first, length))

  private[signal4s] def unsafe(first: Int, length: Int): LagAxis =
    new LagAxis(first, length)

enum CorrelationNormalization:
  case Raw
  case Biased
  case Unbiased
  case Coefficient

/** Cross-correlation result with an explicit lag axis.
  *
  * Convention: \( r_{xy}[\ell] = \sum_n x[n]\, y[n+\ell] \) (raw), then optionally
  * normalized.
  */
final case class Correlation private (
    values: DVec,
    lags: LagAxis
)

object Correlation:
  def apply(values: DVec, lags: LagAxis): Either[SignalError, Correlation] =
    if values.length != lags.length then
      Left(SignalError.LengthMismatch(lags.length, values.length))
    else Right(new Correlation(values, lags))

object Correlate:

  /** Full cross-correlation of `x` with `y` under the documented lag convention. */
  def apply(
      x: DVec,
      y: DVec,
      normalization: CorrelationNormalization = CorrelationNormalization.Raw
  ): Either[SignalError, Correlation] =
    if x.length == 0 || y.length == 0 then Left(SignalError.EmptySignal)
    else
      val n = x.length
      val m = y.length
      // ℓ ranges so that some n has both x[n] and y[n+ℓ] in range:
      // ℓ_min = 1-m, ℓ_max = n-1
      val firstLag = 1 - m
      val lastLag = n - 1
      val outLen = lastLag - firstLag + 1
      val raw = DVecBuilder.zeros(outLen)
      var i = 0
      while i < outLen do
        val lag = firstLag + i
        raw(i) = rawCorrAt(x, y, lag)
        i += 1
      val values = normalize(raw.result(), x, y, normalization)
      Correlation(values, LagAxis.unsafe(firstLag, outLen))

  def auto(
      x: DVec,
      normalization: CorrelationNormalization = CorrelationNormalization.Raw
  ): Either[SignalError, Correlation] =
    apply(x, x, normalization)

  private def rawCorrAt(x: DVec, y: DVec, lag: Int): Double =
    val n = x.length
    val m = y.length
    // n_x in [0, n) and n_x+lag in [0, m)
    val start = math.max(0, -lag)
    val end = math.min(n, m - lag)
    var acc = 0.0
    var t = start
    while t < end do
      acc += x(t) * y(t + lag)
      t += 1
    acc

  private def normalize(
      raw: DVec,
      x: DVec,
      y: DVec,
      normalization: CorrelationNormalization
  ): DVec =
    normalization match
      case CorrelationNormalization.Raw => raw
      case CorrelationNormalization.Biased =>
        scale(raw, 1.0 / x.length.toDouble)
      case CorrelationNormalization.Unbiased =>
        val out = DVecBuilder.zeros(raw.length)
        val n = x.length
        val m = y.length
        val firstLag = 1 - m
        var i = 0
        while i < raw.length do
          val lag = firstLag + i
          val start = math.max(0, -lag)
          val end = math.min(n, m - lag)
          val count = end - start
          out(i) = if count > 0 then raw(i) / count.toDouble else 0.0
          i += 1
        out.result()
      case CorrelationNormalization.Coefficient =>
        val denom = math.sqrt(energy(x) * energy(y))
        if denom == 0.0 then DVec.zeros(raw.length)
        else scale(raw, 1.0 / denom)

  private def scale(v: DVec, alpha: Double): DVec =
    val out = DVecBuilder.zeros(v.length)
    var i = 0
    while i < v.length do
      out(i) = alpha * v(i)
      i += 1
    out.result()

  private def energy(v: DVec): Double =
    var acc = 0.0
    var i = 0
    while i < v.length do
      acc += v(i) * v(i)
      i += 1
    acc
