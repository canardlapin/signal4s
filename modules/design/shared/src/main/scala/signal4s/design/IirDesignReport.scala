package signal4s.design

import signal4s.SampleRate

/** Diagnostics for a designed IIR — warnings, not silent failures. */
final case class IirDesignReport(
    family: String,
    order: Int,
    sampleRate: SampleRate,
    band: FilterBand,
    maxPoleMagnitude: Double,
    stable: Boolean,
    warnings: List[String]
):
  def cutoffSummary: String =
    band match
      case FilterBand.LowPass(c)        => s"lowpass ${c.hertz} Hz"
      case FilterBand.HighPass(c)       => s"highpass ${c.hertz} Hz"
      case FilterBand.BandPass(lo, hi)  => s"bandpass ${lo.hertz}-${hi.hertz} Hz"
      case FilterBand.BandStop(lo, hi)  => s"bandstop ${lo.hertz}-${hi.hertz} Hz"
