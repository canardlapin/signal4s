package signal4s.fft

import gale.linalg.DVec
import signal4s.SignalError

/** STFT result: one complex spectrum column per frame (rfft / onesided layout). */
final case class TimeFrequency private (
    columns: IndexedSeq[ComplexVector],
    times: TimeAxis,
    frequencies: FrequencyAxis,
    nfft: Int,
    onesided: Boolean
):
  def frameCount: Int = columns.length
  def binCount: Int = if columns.isEmpty then 0 else columns.head.length

  /** Magnitude-squared spectrogram (no scaling). */
  def spectrogram: Spectrogram =
    val power = columns.map { col =>
      DVec.tabulate(col.length)(k =>
        val re = col.real(k)
        val im = col.imaginary(k)
        re * re + im * im
      )
    }
    Spectrogram(power, times, frequencies)

object TimeFrequency:
  def apply(
      columns: IndexedSeq[ComplexVector],
      times: TimeAxis,
      frequencies: FrequencyAxis,
      nfft: Int,
      onesided: Boolean
  ): Either[SignalError, TimeFrequency] =
    if times.length != columns.length then
      Left(SignalError.LengthMismatch(times.length, columns.length))
    else if columns.nonEmpty && columns.exists(_.length != frequencies.length) then
      Left(SignalError.LengthMismatch(frequencies.length, columns.head.length))
    else if nfft <= 0 then Left(SignalError.InvalidInputLength(nfft))
    else Right(new TimeFrequency(columns, times, frequencies, nfft, onesided))

/** Derived power view of an STFT — not an independent transform. */
final case class Spectrogram(
    power: IndexedSeq[DVec],
    times: TimeAxis,
    frequencies: FrequencyAxis
):
  def frameCount: Int = power.length
