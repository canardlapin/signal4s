package signal4s.fft

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*

/** Welch PSD / power-spectrum estimator sharing [[FramePlan]]. */
final class WelchPlan private (
    val window: Window,
    val frames: FramePlan,
    val nfft: Int,
    val sampleRate: SampleRate,
    val detrend: Detrend,
    val sides: SpectralSides,
    val scaling: SpectralScaling,
    val average: AverageMethod,
    private val fft: RealFftPlan,
    private val sumW: Double,
    private val sumW2: Double
):
  def estimate(signal: DVec): Either[SignalError, WelchResult] =
    sides match
      case SpectralSides.Twosided =>
        Left(SignalError.NumericalFailure("WelchPlan", "twosided not implemented yet"))
      case SpectralSides.Onesided =>
        // SciPy welch: only segments with start + nperseg <= n (no trailing pad).
        val nFull = welchSegmentCount(signal.length)
        if nFull <= 0 then
          Left(SignalError.NumericalFailure("WelchPlan", "signal shorter than one segment"))
        else
          val binCount = nfft / 2 + 1
          val acc = Array.fill(binCount)(0.0)
          val accSq = if average == AverageMethod.Median then Array.ofDim[Double](nFull * binCount) else null
          val workspace = fft.newWorkspace()
          var s = 0
          var err: Option[SignalError] = None
          while s < nFull && err.isEmpty do
            extractWelchFrame(signal, s) match
              case Left(e) => err = Some(e)
              case Right(raw) =>
                val detrended = DetrendOps(raw, detrend)
                val windowed = applyWindow(detrended)
                fft.forwardInto(pad(windowed), workspace) match
                  case Left(e) => err = Some(e)
                  case Right(_) =>
                    var k = 0
                    while k < binCount do
                      val re = workspace.re(k)
                      val im = workspace.im(k)
                      val p = re * re + im * im
                      average match
                        case AverageMethod.Mean =>
                          acc(k) += p
                        case AverageMethod.Median =>
                          accSq(s * binCount + k) = p
                      k += 1
            s += 1
          err match
            case Some(e) => Left(e)
            case None =>
              val power = DVecBuilder.zeros(binCount)
              average match
                case AverageMethod.Mean =>
                  var k = 0
                  while k < binCount do
                    power(k) = scaleBin(acc(k) / nFull.toDouble, k, binCount)
                    k += 1
                case AverageMethod.Median =>
                  var k = 0
                  while k < binCount do
                    val vals = Array.ofDim[Double](nFull)
                    var s2 = 0
                    while s2 < nFull do
                      vals(s2) = accSq(s2 * binCount + k)
                      s2 += 1
                    java.util.Arrays.sort(vals)
                    val med =
                      if nFull % 2 == 1 then vals(nFull / 2)
                      else 0.5 * (vals(nFull / 2 - 1) + vals(nFull / 2))
                    power(k) = scaleBin(med, k, binCount)
                    k += 1
              FrequencyAxis.realFft(nfft, sampleRate).map { freqs =>
                WelchResult(
                  frequencies = freqs,
                  power = power.result(),
                  scaling = scaling,
                  sides = sides,
                  segmentCount = nFull,
                  windowPower = sumW2,
                  windowCoherentPower = sumW * sumW,
                  detrend = detrend,
                  average = average,
                  degreesOfFreedom = None,
                  diagnostics = List(
                    s"hop=${frames.hop}",
                    s"nfft=$nfft",
                    "dof=unspecified (not fabricated)"
                  )
                )
              }

  private def welchSegmentCount(signalLength: Int): Int =
    if signalLength < frames.frameLength then 0
    else 1 + (signalLength - frames.frameLength) / frames.hop

  private def extractWelchFrame(signal: DVec, segment: Int): Either[SignalError, DVec] =
    val start = segment * frames.hop
    if start + frames.frameLength > signal.length then
      Left(SignalError.NumericalFailure("WelchPlan", s"segment $segment out of range"))
    else
      val out = DVecBuilder.zeros(frames.frameLength)
      var i = 0
      while i < frames.frameLength do
        out(i) = signal(start + i)
        i += 1
      Right(out.result())

  private def applyWindow(frame: DVec): DVec =
    DVec.tabulate(frame.length)(i => frame(i) * window.taps(i))

  private def pad(frame: DVec): DVec =
    if frame.length == nfft then frame
    else
      val out = DVecBuilder.zeros(nfft)
      var i = 0
      while i < frame.length do
        out(i) = frame(i)
        i += 1
      out.result()

  /** SciPy onesided scaling for density / spectrum. */
  private def scaleBin(meanPower: Double, k: Int, binCount: Int): Double =
    val fs = sampleRate.hertz
    val onesidedFactor =
      if k == 0 || (nfft % 2 == 0 && k == binCount - 1) then 1.0
      else 2.0
    scaling match
      case SpectralScaling.Density =>
        onesidedFactor * meanPower / (fs * sumW2)
      case SpectralScaling.Spectrum =>
        onesidedFactor * meanPower / (sumW * sumW)

object WelchPlan:
  def apply(
      window: Window,
      hop: Int,
      sampleRate: SampleRate,
      nfft: Int,
      detrend: Detrend = Detrend.Mean,
      sides: SpectralSides = SpectralSides.Onesided,
      scaling: SpectralScaling = SpectralScaling.Density,
      average: AverageMethod = AverageMethod.Mean
  ): Either[SignalError, WelchPlan] =
    if nfft < window.length then
      Left(SignalError.NumericalFailure("WelchPlan", s"nfft=$nfft < window ${window.length}"))
    else
      FramePlan(window.length, hop, FrameAlignment.FromStart, Boundary.Zero).flatMap { frames =>
        RealFftPlan(nfft, FftNormalization.Backward, sampleRate).map { fft =>
          var sumW = 0.0
          var sumW2 = 0.0
          var i = 0
          while i < window.length do
            val w = window.taps(i)
            sumW += w
            sumW2 += w * w
            i += 1
          new WelchPlan(
            window,
            frames,
            nfft,
            sampleRate,
            detrend,
            sides,
            scaling,
            average,
            fft,
            sumW,
            sumW2
          )
        }
      }

  def hann(
      nperseg: Int,
      sampleRate: SampleRate,
      noverlap: Option[Int] = None,
      nfft: Option[Int] = None,
      detrend: Detrend = Detrend.Mean,
      scaling: SpectralScaling = SpectralScaling.Density
  ): Either[SignalError, WelchPlan] =
    val overlap = noverlap.getOrElse(nperseg / 2)
    val hop = nperseg - overlap
    val n = nfft.getOrElse(nperseg)
    Window
      .fromSpec(WindowSpec.Hann(nperseg, WindowConvention.Periodic))
      .flatMap(w => apply(w, hop, sampleRate, n, detrend = detrend, scaling = scaling))
