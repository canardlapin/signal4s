package signal4s.fft

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*

/** Analysis/synthesis STFT plan with invertibility validated at construction.
  *
  * Frame placement matches SciPy `ShortTimeFFT` for even windows with
  * [[FrameAlignment.Centered]] (first frame starts at `-frameLength/2`).
  */
final class StftPlan private (
    val analysis: Window,
    val synthesis: Window,
    val frames: FramePlan,
    val nfft: Int,
    val sampleRate: SampleRate,
    val dual: Option[CanonicalDual],
    private val fft: RealFftPlan
):
  def hop: Int = frames.hop
  def frameLength: Int = frames.frameLength

  def analyze(signal: DVec): Either[SignalError, TimeFrequency] =
    val nFrames = frames.frameCount(signal.length)
    if nFrames == 0 then
      Left(SignalError.NumericalFailure("StftPlan.analyze", "no frames for signal length"))
    else
      val cols = Vector.newBuilder[ComplexVector]
      val workspace = fft.newWorkspace()
      var p = 0
      var err: Option[SignalError] = None
      while p < nFrames && err.isEmpty do
        frames.extract(signal, p, analysis) match
          case Left(e) => err = Some(e)
          case Right(block) =>
            fft.forwardInto(padToNfft(block), workspace) match
              case Left(e) => err = Some(e)
              case Right(_) =>
                fft.spectrumFrom(workspace) match
                  case Left(e)     => err = Some(e)
                  case Right(spec) => cols += spec.bins
        p += 1
      err match
        case Some(e) => Left(e)
        case None =>
          for
            times <- TimeAxis.fromFrames(nFrames, hop, sampleRate)
            freqs <- FrequencyAxis.realFft(nfft, sampleRate)
            tf <- TimeFrequency(cols.result(), times, freqs, nfft, onesided = true)
          yield tf

  def synthesize(tf: TimeFrequency, outputLength: Int): Either[SignalError, DVec] =
    if tf.nfft != nfft then Left(SignalError.LengthMismatch(nfft, tf.nfft))
    else if !tf.onesided then
      Left(SignalError.NumericalFailure("StftPlan.synthesize", "twosided STFT not supported yet"))
    else if outputLength <= 0 then Left(SignalError.InvalidInputLength(outputLength))
    else
      val synWin = dual.map(_.taps).getOrElse(synthesis.taps)
      val dest = Array.fill(outputLength)(0.0)
      val workspace = fft.newWorkspace()
      var p = 0
      var err: Option[SignalError] = None
      while p < tf.frameCount && err.isEmpty do
        val col = tf.columns(p)
        RealSpectrum(col, tf.frequencies, nfft, FftNormalization.Backward) match
          case Left(e) => err = Some(e)
          case Right(spec) =>
            fft.inverseInto(spec, workspace) match
              case Left(e) => err = Some(e)
              case Right(_) =>
                val framed = DVecBuilder.zeros(frameLength)
                var i = 0
                while i < frameLength do
                  framed(i) = workspace.re(i) * synWin(i)
                  i += 1
                frames.overlapAdd(dest, framed.result(), p)
        p += 1
      err match
        case Some(e) => Left(e)
        case None =>
          val out = DVecBuilder.zeros(outputLength)
          var i = 0
          while i < outputLength do
            out(i) = dest(i)
            i += 1
          Right(out.result())

  private def padToNfft(block: DVec): DVec =
    if block.length == nfft then block
    else
      val out = DVecBuilder.zeros(nfft)
      var i = 0
      while i < block.length do
        out(i) = block(i)
        i += 1
      out.result()

object StftPlan:
  def apply(
      analysis: Window,
      hop: Int,
      sampleRate: SampleRate,
      nfft: Int,
      alignment: FrameAlignment = FrameAlignment.Centered,
      boundary: Boundary = Boundary.Zero,
      synthesis: Option[Window] = None
  ): Either[SignalError, StftPlan] =
    if nfft < analysis.length then
      Left(
        SignalError.NumericalFailure(
          "StftPlan",
          s"nfft=$nfft must be >= frame length ${analysis.length}"
        )
      )
    else
      FramePlan(analysis.length, hop, alignment, boundary).flatMap { frames =>
        RealFftPlan(nfft, FftNormalization.Backward, sampleRate).flatMap { fft =>
          synthesis match
            case Some(syn) =>
              if syn.length != analysis.length then
                Left(SignalError.LengthMismatch(analysis.length, syn.length))
              else
                Right(new StftPlan(analysis, syn, frames, nfft, sampleRate, None, fft))
            case None =>
              CanonicalDual.fromAnalysis(analysis, hop).flatMap { d =>
                Window.fromTaps(d.taps, analysis.convention).map { synWin =>
                  new StftPlan(analysis, synWin, frames, nfft, sampleRate, Some(d), fft)
                }
              }
        }
      }

  /** Convenience: Hann periodic window, 50% hop, centered frames. */
  def hannPeriodic(
      frameLength: Int,
      sampleRate: SampleRate,
      hop: Option[Int] = None,
      nfft: Option[Int] = None
  ): Either[SignalError, StftPlan] =
    val h = hop.getOrElse(frameLength / 2)
    val n = nfft.getOrElse(frameLength)
    Window
      .fromSpec(WindowSpec.Hann(frameLength, WindowConvention.Periodic))
      .flatMap(win => apply(win, h, sampleRate, n))
