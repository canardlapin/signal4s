package signal4s

/** Sampling frequency in hertz. Must be finite and strictly positive. */
final case class SampleRate private (hertz: Double):
  def nyquist: Frequency = Frequency.unsafe(hertz / 2.0)

  def cyclesPerSample(frequency: Frequency): Either[SignalError, CyclesPerSample] =
    frequency.cyclesPerSample(at = this)

  def radiansPerSample(frequency: Frequency): Either[SignalError, RadiansPerSample] =
    frequency.radiansPerSample(at = this)

object SampleRate:
  def hertz(value: Double): Either[SignalError, SampleRate] =
    if value.isFinite && value > 0.0 then Right(SampleRate(value))
    else Left(SignalError.InvalidSampleRate(value))

  private[signal4s] def unsafe(value: Double): SampleRate = SampleRate(value)

/** Physical frequency in hertz. Must be finite and non-negative. */
final case class Frequency private (hertz: Double):
  def cyclesPerSample(at: SampleRate): Either[SignalError, CyclesPerSample] =
    if hertz > at.nyquist.hertz then
      Left(SignalError.FrequencyAboveNyquist(this, at))
    else Right(CyclesPerSample.unsafe(hertz / at.hertz))

  def radiansPerSample(at: SampleRate): Either[SignalError, RadiansPerSample] =
    cyclesPerSample(at).map(c => RadiansPerSample.unsafe(2.0 * math.Pi * c.value))

object Frequency:
  def hertz(value: Double): Either[SignalError, Frequency] =
    if value.isFinite && value >= 0.0 then Right(Frequency(value))
    else Left(SignalError.InvalidFrequency(value))

  private[signal4s] def unsafe(value: Double): Frequency = Frequency(value)

/** Normalized frequency in cycles per sample (1.0 = sampling rate). */
final case class CyclesPerSample private (value: Double):
  def toFrequency(at: SampleRate): Frequency =
    Frequency.unsafe(value * at.hertz)

  def toRadiansPerSample: RadiansPerSample =
    RadiansPerSample.unsafe(2.0 * math.Pi * value)

object CyclesPerSample:
  def of(value: Double): Either[SignalError, CyclesPerSample] =
    if value.isFinite && value >= 0.0 && value <= 0.5 then Right(unsafe(value))
    else Left(SignalError.InvalidFrequency(value))

  private[signal4s] def unsafe(value: Double): CyclesPerSample =
    new CyclesPerSample(value)

/** Normalized frequency in radians per sample (π = Nyquist). */
final case class RadiansPerSample private (value: Double):
  def toCyclesPerSample: CyclesPerSample =
    CyclesPerSample.unsafe(value / (2.0 * math.Pi))

  def toFrequency(at: SampleRate): Frequency =
    toCyclesPerSample.toFrequency(at)

object RadiansPerSample:
  def of(value: Double): Either[SignalError, RadiansPerSample] =
    if value.isFinite && value >= 0.0 && value <= math.Pi then Right(unsafe(value))
    else Left(SignalError.InvalidFrequency(value))

  private[signal4s] def unsafe(value: Double): RadiansPerSample =
    new RadiansPerSample(value)

/** Time offset or duration in seconds. Must be finite. */
final case class Seconds private (value: Double):
  def +(other: Seconds): Seconds = Seconds.unsafe(value + other.value)
  def -(other: Seconds): Seconds = Seconds.unsafe(value - other.value)

object Seconds:
  val Zero: Seconds = Seconds.unsafe(0.0)

  def of(value: Double): Either[SignalError, Seconds] =
    if value.isFinite then Right(Seconds.unsafe(value))
    else Left(SignalError.InvalidDuration(value))

  private[signal4s] def unsafe(value: Double): Seconds = Seconds(value)
