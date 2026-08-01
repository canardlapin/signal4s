package signal4s.fft

import signal4s.{SampleRate, Seconds, SignalError}

/** Regular time grid for STFT frames (seconds). */
final case class TimeAxis private (
    first: Seconds,
    step: Seconds,
    length: Int
):
  def timeAt(index: Int): Seconds =
    require(index >= 0 && index < length, s"time index $index out of range")
    Seconds.unsafe(first.value + index.toDouble * step.value)

object TimeAxis:
  def apply(first: Seconds, step: Seconds, length: Int): Either[SignalError, TimeAxis] =
    if length < 0 then Left(SignalError.InvalidInputLength(length))
    else if !step.value.isFinite then Left(SignalError.InvalidDuration(step.value))
    else Right(new TimeAxis(first, step, length))

  def fromFrames(
      frameCount: Int,
      hop: Int,
      sampleRate: SampleRate,
      firstFrameCenterSample: Double = 0.0
  ): Either[SignalError, TimeAxis] =
    for
      first <- Seconds.of(firstFrameCenterSample / sampleRate.hertz)
      step <- Seconds.of(hop.toDouble / sampleRate.hertz)
      axis <- apply(first, step, frameCount)
    yield axis
