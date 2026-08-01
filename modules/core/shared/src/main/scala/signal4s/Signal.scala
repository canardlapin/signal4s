package signal4s

import gale.linalg.DVec

/** Finite observation of a discrete-time signal on a regular sampling lattice. */
final case class Signal private (
    samples: DVec,
    sampling: Sampling
):
  def length: Int = samples.length
  def sampleRate: SampleRate = sampling.rate
  def start: Seconds = sampling.start

  def timeAt(index: Int): Seconds = sampling.timeAt(index)

object Signal:
  def apply(
      samples: DVec,
      sampling: Sampling
  ): Either[SignalError, Signal] =
    if samples.length <= 0 then Left(SignalError.EmptySignal)
    else Right(new Signal(samples, sampling))
