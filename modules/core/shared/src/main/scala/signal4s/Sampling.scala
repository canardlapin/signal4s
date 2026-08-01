package signal4s

/** Regular sampling lattice: rate and absolute time of the first sample. */
final case class Sampling private (rate: SampleRate, start: Seconds):
  def samplePeriod: Seconds = Seconds.unsafe(1.0 / rate.hertz)

  def timeAt(index: Int): Seconds =
    Seconds.unsafe(start.value + index.toDouble / rate.hertz)

  def shiftBySamples(samples: Int): Sampling =
    Sampling(
      rate,
      Seconds.unsafe(start.value + samples.toDouble / rate.hertz)
    )

object Sampling:
  def apply(rate: SampleRate, start: Seconds = Seconds.Zero): Sampling =
    new Sampling(rate, start)
