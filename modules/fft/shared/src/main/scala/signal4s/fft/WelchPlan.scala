package signal4s.fft

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*
import gale.numeric.ExactSum

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
          val acc = Array.fill(binCount)(ExactSum.zero())
          val normalizer = ExactSum.zero()
          val _ = normalizer.add(nFull.toDouble)
          if average == AverageMethod.Median && BigInt(nFull) * binCount > Int.MaxValue then
            return Left(SignalError.NumericalFailure("WelchPlan", "median workspace exceeds Int capacity"))
          val accSq = if average == AverageMethod.Median then Array.ofDim[Double](nFull * binCount) else null
          val workspace = fft.newWorkspace()
          var s = 0
          var err: Option[SignalError] = None
          while s < nFull && err.isEmpty do
            extractWelchFrame(signal, s) match
              case Left(e) => err = Some(e)
              case Right(raw) =>
                SpectralFrame.prepare(raw, window, detrend).flatMap { case (windowed, scale) =>
                  fft.forwardInto(pad(windowed), workspace).map(_ => scale)
                } match
                  case Left(e) => err = Some(e)
                  case Right(scale) =>
                    var k = 0
                    while k < binCount && err.isEmpty do
                      val re = workspace.re(k)
                      val im = workspace.im(k)
                      val factor = if k == 0 || (nfft % 2 == 0 && k == binCount - 1) then 1.0 else 2.0
                      val divisor = if scaling == SpectralScaling.Density then sampleRate.hertz * sumW2 else sumW * sumW
                      SpectralFrame.power(re, im, scale, divisor, factor) match
                        case Left(e) => err = Some(e)
                        case Right(p) => average match
                          case AverageMethod.Mean =>
                            acc(k).add(p) match
                              case Left(e) => err = Some(SignalError.NumericalFailure("WelchPlan", e.message))
                              case Right(_) => ()
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
                    acc(k).ratio(normalizer) match
                      case Left(e) => return Left(SignalError.NumericalFailure("WelchPlan", e.message))
                      case Right(value) => power(k) = value
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
                      else
                        val middle = ExactSum.zero()
                        val two = ExactSum.zero()
                        val _ = middle.add(vals(nFull / 2 - 1))
                        val _ = middle.add(vals(nFull / 2))
                        val _ = two.add(2.0)
                        middle.ratio(two) match
                          case Left(e) => return Left(SignalError.NumericalFailure("WelchPlan", e.message))
                          case Right(value) => value
                    power(k) = med
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

  private def pad(frame: DVec): DVec =
    if frame.length == nfft then frame
    else
      val out = DVecBuilder.zeros(nfft)
      var i = 0
      while i < frame.length do
        out(i) = frame(i)
        i += 1
      out.result()

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
        RealFftPlan(nfft, FftNormalization.Backward, sampleRate).flatMap { fft =>
          var sumW = 0.0
          var sumW2 = 0.0
          var i = 0
          while i < window.length do
            val w = window.taps(i)
            sumW += w
            sumW2 += w * w
            i += 1
          val normalization = scaling match
            case SpectralScaling.Density => sampleRate.hertz * sumW2
            case SpectralScaling.Spectrum => sumW * sumW
          if !normalization.isFinite || normalization <= 0.0 then
            // A plan must not publish successful NaN/Infinity from an undefined divisor.
            Left(SignalError.NumericalFailure("WelchPlan", "spectral normalization must be finite and positive"))
          else Right(new WelchPlan(
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
          ))
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
