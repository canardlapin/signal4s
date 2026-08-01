package signal4s.design

import signal4s.*

/** Digital Butterworth design via analog prototype → prewarp → bilinear → SOS. */
object Butterworth:

  def lowPass(
      order: Int,
      cutoff: Frequency,
      sampleRate: SampleRate
  ): Either[SignalError, DesignedIir] =
    design(order, FilterBand.LowPass(cutoff), sampleRate)

  def highPass(
      order: Int,
      cutoff: Frequency,
      sampleRate: SampleRate
  ): Either[SignalError, DesignedIir] =
    design(order, FilterBand.HighPass(cutoff), sampleRate)

  def bandPass(
      order: Int,
      low: Frequency,
      high: Frequency,
      sampleRate: SampleRate
  ): Either[SignalError, DesignedIir] =
    design(order, FilterBand.BandPass(low, high), sampleRate)

  def design(
      order: Int,
      band: FilterBand,
      sampleRate: SampleRate
  ): Either[SignalError, DesignedIir] =
    if order <= 0 then Left(SignalError.InvalidInputLength(order))
    else
      for
        _ <- band.validate(sampleRate)
        proto <- AnalogPrototype.butterworth(order)
        analog <- frequencyTransform(proto, band, sampleRate)
        digital <- Bilinear.transform(analog, sampleRate)
        sos <- Sos.fromZerosPolesGain(digital)
        report <- Right(buildReport(order, band, sampleRate, digital))
      yield DesignedIir(digital, sos, report)

  /** Prewarped analog edge frequencies in rad/s: \( \Omega = 2 f_s \tan(\pi f / f_s) \). */
  private def prewarp(f: Frequency, sampleRate: SampleRate): Double =
    2.0 * sampleRate.hertz * math.tan(math.Pi * f.hertz / sampleRate.hertz)

  private def frequencyTransform(
      proto: ZerosPolesGain,
      band: FilterBand,
      sampleRate: SampleRate
  ): Either[SignalError, ZerosPolesGain] =
    band match
      case FilterBand.LowPass(c) =>
        AnalogPrototype.lp2lp(proto, prewarp(c, sampleRate))
      case FilterBand.HighPass(c) =>
        AnalogPrototype.lp2hp(proto, prewarp(c, sampleRate))
      case FilterBand.BandPass(lo, hi) =>
        val w0 = math.sqrt(prewarp(lo, sampleRate) * prewarp(hi, sampleRate))
        val bw = prewarp(hi, sampleRate) - prewarp(lo, sampleRate)
        AnalogPrototype.lp2bp(proto, w0, bw)
      case FilterBand.BandStop(_, _) =>
        Left(
          SignalError.NumericalFailure(
            "Butterworth",
            "BandStop not implemented in this slice (lp2bs deferred)"
          )
        )

  private def buildReport(
      order: Int,
      band: FilterBand,
      sampleRate: SampleRate,
      zpk: ZerosPolesGain
  ): IirDesignReport =
    var maxMag = 0.0
    var i = 0
    while i < zpk.poles.length do
      maxMag = math.max(maxMag, zpk.poles(i).abs)
      i += 1
    val stable = maxMag < 1.0 - 1e-12
    val warnings = List.newBuilder[String]
    if maxMag > 0.98 then
      warnings += f"max |pole|=$maxMag%.6f is close to the unit circle"
    if order >= 8 then
      warnings += "high-order IIR: prefer SOS execution over expanded TF"
    if !stable then warnings += "design is unstable (|pole| >= 1)"
    IirDesignReport(
      family = "Butterworth",
      order = order,
      sampleRate = sampleRate,
      band = band,
      maxPoleMagnitude = maxMag,
      stable = stable,
      warnings = warnings.result()
    )
