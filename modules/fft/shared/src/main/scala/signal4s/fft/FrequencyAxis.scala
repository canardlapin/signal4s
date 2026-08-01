package signal4s.fft

import signal4s.{Frequency, SampleRate, SignalError}

/** Regular frequency grid for spectra. */
final case class FrequencyAxis private (
    first: Frequency,
    step: Frequency,
    length: Int
):
  def frequencyAt(index: Int): Frequency =
    require(index >= 0 && index < length, s"frequency index $index out of range")
    Frequency.unsafe(first.hertz + index.toDouble * step.hertz)

object FrequencyAxis:
  def apply(
      first: Frequency,
      step: Frequency,
      length: Int
  ): Either[SignalError, FrequencyAxis] =
    if length < 0 then Left(SignalError.InvalidInputLength(length))
    else if !step.hertz.isFinite then Left(SignalError.InvalidFrequency(step.hertz))
    else Right(new FrequencyAxis(first, step, length))

  /** Hertz axis for a real FFT of `sourceLength` at `sampleRate` (rfft bins). */
  def realFft(
      sourceLength: Int,
      sampleRate: SampleRate
  ): Either[SignalError, FrequencyAxis] =
    if sourceLength <= 0 then Left(SignalError.InvalidInputLength(sourceLength))
    else
      val bins = sourceLength / 2 + 1
      val stepHz = sampleRate.hertz / sourceLength.toDouble
      for
        first <- Frequency.hertz(0.0)
        step <- Frequency.hertz(stepHz)
        axis <- apply(first, step, bins)
      yield axis

  /** Normalized cycles/sample axis stored as Frequency with fs = 1 Hz. */
  def normalizedRealFft(sourceLength: Int): Either[SignalError, FrequencyAxis] =
    realFft(sourceLength, SampleRate.unsafe(1.0))
