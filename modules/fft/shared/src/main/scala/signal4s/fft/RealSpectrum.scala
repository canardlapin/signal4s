package signal4s.fft

import signal4s.SignalError

/** Real-input FFT result: `bins.length == sourceLength/2 + 1`. */
final case class RealSpectrum private (
    bins: ComplexVector,
    frequencies: FrequencyAxis,
    sourceLength: Int,
    normalization: FftNormalization
):
  def length: Int = bins.length

object RealSpectrum:
  def apply(
      bins: ComplexVector,
      frequencies: FrequencyAxis,
      sourceLength: Int,
      normalization: FftNormalization
  ): Either[SignalError, RealSpectrum] =
    val expected = sourceLength / 2 + 1
    if sourceLength <= 0 then Left(SignalError.InvalidInputLength(sourceLength))
    else if bins.length != expected then
      Left(SignalError.LengthMismatch(expected, bins.length))
    else if frequencies.length != expected then
      Left(SignalError.LengthMismatch(expected, frequencies.length))
    else Right(new RealSpectrum(bins, frequencies, sourceLength, normalization))
