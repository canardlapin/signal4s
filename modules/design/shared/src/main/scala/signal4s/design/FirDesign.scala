package signal4s.design

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*
import signal4s.filter.Fir
import signal4s.fft.{Window, WindowConvention, WindowSpec}

/** Windowed-sinc FIR design (SciPy `firwin` style) with physical frequencies. */
object FirDesign:

  def lowPass(
      numTaps: Int,
      cutoff: Frequency,
      sampleRate: SampleRate,
      window: WindowSpec = WindowSpec.Hamming(1, WindowConvention.Symmetric)
  ): Either[SignalError, DesignedFir] =
    design(numTaps, FilterBand.LowPass(cutoff), sampleRate, window)

  def highPass(
      numTaps: Int,
      cutoff: Frequency,
      sampleRate: SampleRate,
      window: WindowSpec = WindowSpec.Hamming(1, WindowConvention.Symmetric)
  ): Either[SignalError, DesignedFir] =
    design(numTaps, FilterBand.HighPass(cutoff), sampleRate, window)

  def bandPass(
      numTaps: Int,
      low: Frequency,
      high: Frequency,
      sampleRate: SampleRate,
      window: WindowSpec = WindowSpec.Hamming(1, WindowConvention.Symmetric)
  ): Either[SignalError, DesignedFir] =
    design(numTaps, FilterBand.BandPass(low, high), sampleRate, window)

  def design(
      numTaps: Int,
      band: FilterBand,
      sampleRate: SampleRate,
      windowSpec: WindowSpec
  ): Either[SignalError, DesignedFir] =
    if numTaps <= 0 then Left(SignalError.InvalidInputLength(numTaps))
    else if numTaps % 2 == 0 && (band.isInstanceOf[FilterBand.HighPass] ||
        band.isInstanceOf[FilterBand.BandPass])
    then
      Left(
        SignalError.NumericalFailure(
          "FirDesign",
          "highpass/bandpass firwin-style design requires an odd number of taps"
        )
      )
    else
      for
        _ <- band.validate(sampleRate)
        win <- Window.fromSpec(withLength(windowSpec, numTaps))
        taps <- idealWindowed(numTaps, band, sampleRate, win)
        fir <- Fir.causal(taps)
      yield DesignedFir(fir, win, band, numTaps)

  /** Kaiser order estimate (`scipy.signal.kaiserord`).
    *
    * @param rippleDb stopband attenuation in dB (positive)
    * @param transitionWidthHz width of transition band in hertz
    */
  def kaiserOrder(
      rippleDb: Double,
      transitionWidthHz: Double,
      sampleRate: SampleRate
  ): Either[SignalError, (Int, Double)] =
    if !(rippleDb > 0.0 && rippleDb.isFinite) then
      Left(SignalError.NumericalFailure("kaiserOrder", s"rippleDb must be > 0, got $rippleDb"))
    else if !(transitionWidthHz > 0.0 && transitionWidthHz.isFinite) then
      Left(SignalError.NumericalFailure("kaiserOrder", "transition width must be positive"))
    else if transitionWidthHz >= sampleRate.nyquist.hertz then
      Left(SignalError.FrequencyAboveNyquist(Frequency.unsafe(transitionWidthHz), sampleRate))
    else
      val width = transitionWidthHz / sampleRate.nyquist.hertz
      val a = rippleDb
      val beta =
        if a > 50 then 0.1102 * (a - 8.7)
        else if a > 21 then 0.5842 * math.pow(a - 21, 0.4) + 0.07886 * (a - 21)
        else 0.0
      val numtaps = ((a - 8.0) / (2.285 * math.Pi * width)).ceil.toInt + 1
      Right((math.max(1, numtaps), beta))

  private def withLength(spec: WindowSpec, n: Int): WindowSpec =
    spec match
      case WindowSpec.Rectangular(_, c)  => WindowSpec.Rectangular(n, c)
      case WindowSpec.Hann(_, c)         => WindowSpec.Hann(n, c)
      case WindowSpec.Hamming(_, c)      => WindowSpec.Hamming(n, c)
      case WindowSpec.Blackman(_, c)     => WindowSpec.Blackman(n, c)
      case WindowSpec.Kaiser(_, beta, c) => WindowSpec.Kaiser(n, beta, c)

  private def idealWindowed(
      numTaps: Int,
      band: FilterBand,
      sampleRate: SampleRate,
      window: Window
  ): Either[SignalError, DVec] =
    val m = (numTaps - 1).toDouble / 2.0
    val fs = sampleRate.hertz
    val ideal = DVecBuilder.zeros(numTaps)
    var n = 0
    while n < numTaps do
      val k = n.toDouble - m
      ideal(n) = band match
        case FilterBand.LowPass(c) =>
          2.0 * c.hertz / fs * sinc(2.0 * c.hertz / fs * k)
        case FilterBand.HighPass(c) =>
          val lp = 2.0 * c.hertz / fs * sinc(2.0 * c.hertz / fs * k)
          if math.abs(k) < 1e-16 then 1.0 - lp else -lp
        case FilterBand.BandPass(lo, hi) =>
          val hHi = 2.0 * hi.hertz / fs * sinc(2.0 * hi.hertz / fs * k)
          val hLo = 2.0 * lo.hertz / fs * sinc(2.0 * lo.hertz / fs * k)
          hHi - hLo
        case FilterBand.BandStop(lo, hi) =>
          val hHi = 2.0 * hi.hertz / fs * sinc(2.0 * hi.hertz / fs * k)
          val hLo = 2.0 * lo.hertz / fs * sinc(2.0 * lo.hertz / fs * k)
          (if math.abs(k) < 1e-16 then 1.0 else 0.0) - (hHi - hLo)
      n += 1
    val out = DVecBuilder.zeros(numTaps)
    n = 0
    while n < numTaps do
      out(n) = ideal(n) * window.taps(n)
      n += 1
    Right(scalePassband(out.result(), band, sampleRate))

  private def scalePassband(h: DVec, band: FilterBand, sampleRate: SampleRate): DVec =
    val scale =
      band match
        case FilterBand.LowPass(_) | FilterBand.BandStop(_, _) =>
          sum(h)
        case FilterBand.HighPass(_) =>
          alternatingSum(h)
        case FilterBand.BandPass(lo, hi) =>
          val f0 = 0.5 * (lo.hertz + hi.hertz)
          val omega = 2.0 * math.Pi * f0 / sampleRate.hertz
          dtftReal(h, omega)
    if math.abs(scale) < 1e-30 then h
    else DVec.tabulate(h.length)(i => h(i) / scale)

  private def sum(h: DVec): Double =
    var s = 0.0
    var i = 0
    while i < h.length do
      s += h(i)
      i += 1
    s

  private def alternatingSum(h: DVec): Double =
    var s = 0.0
    var i = 0
    while i < h.length do
      s += (if i % 2 == 0 then h(i) else -h(i))
      i += 1
    s

  private def dtftReal(h: DVec, omega: Double): Double =
    var re = 0.0
    var i = 0
    while i < h.length do
      re += h(i) * math.cos(omega * i)
      i += 1
    re

  private def sinc(x: Double): Double =
    if math.abs(x) < 1e-16 then 1.0
    else math.sin(math.Pi * x) / (math.Pi * x)
