package signal4s.design

import signal4s.{Frequency, SampleRate, SignalError}

/** Band specification with physical frequencies (never bare cutoffs). */
enum FilterBand:
  case LowPass(cutoff: Frequency)
  case HighPass(cutoff: Frequency)
  case BandPass(low: Frequency, high: Frequency)
  case BandStop(low: Frequency, high: Frequency)

  def validate(sampleRate: SampleRate): Either[SignalError, Unit] =
    this match
      case LowPass(c) =>
        requireBelowNyquist(c, sampleRate)
      case HighPass(c) =>
        requireBelowNyquist(c, sampleRate)
      case BandPass(lo, hi) =>
        for
          _ <- requireBelowNyquist(lo, sampleRate)
          _ <- requireBelowNyquist(hi, sampleRate)
          _ <-
            if lo.hertz >= hi.hertz then
              Left(
                SignalError.NumericalFailure(
                  "FilterBand",
                  s"band edges require low < high, got ${lo.hertz} >= ${hi.hertz}"
                )
              )
            else Right(())
        yield ()
      case BandStop(lo, hi) =>
        for
          _ <- requireBelowNyquist(lo, sampleRate)
          _ <- requireBelowNyquist(hi, sampleRate)
          _ <-
            if lo.hertz >= hi.hertz then
              Left(
                SignalError.NumericalFailure(
                  "FilterBand",
                  s"band edges require low < high, got ${lo.hertz} >= ${hi.hertz}"
                )
              )
            else Right(())
        yield ()

  private def requireBelowNyquist(
      f: Frequency,
      sampleRate: SampleRate
  ): Either[SignalError, Unit] =
    if f.hertz > sampleRate.nyquist.hertz then
      Left(SignalError.FrequencyAboveNyquist(f, sampleRate))
    else if f.hertz <= 0.0 then
      Left(SignalError.InvalidFrequency(f.hertz))
    else Right(())
