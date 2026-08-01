package signal4s.fft

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*
import signal4s.internal.SampleAccess

/** How the first analysis frame is placed relative to sample index 0. */
enum FrameAlignment:
  /** First frame starts at sample 0 (Welch / periodogram segments). */
  case FromStart
  /** First frame starts at `-(frameLength/2)` — SciPy [[ShortTimeFFT]] for even windows. */
  case Centered

/** Shared framing for STFT, Welch, and other block spectral tools. */
final case class FramePlan private (
    frameLength: Int,
    hop: Int,
    alignment: FrameAlignment,
    boundary: Boundary
):
  def firstFrameStart: Int =
    alignment match
      case FrameAlignment.FromStart => 0
      case FrameAlignment.Centered  => -(frameLength / 2)

  /** Number of frames whose start is `< signalLength` (ShortTimeFFT `p_max` style). */
  def frameCount(signalLength: Int): Int =
    if signalLength <= 0 then 0
    else
      val start0 = firstFrameStart
      // last frame index p with start0 + p*hop < signalLength
      val raw = signalLength - start0
      if raw <= 0 then 0
      else (raw + hop - 1) / hop

  def frameStart(frameIndex: Int): Int =
    firstFrameStart + frameIndex * hop

  /** Extract one framed, windowed block (length = frameLength). */
  def extract(
      signal: DVec,
      frameIndex: Int,
      window: Window
  ): Either[SignalError, DVec] =
    if window.length != frameLength then
      Left(SignalError.LengthMismatch(frameLength, window.length))
    else if frameIndex < 0 || frameIndex >= frameCount(signal.length) then
      Left(
        SignalError.NumericalFailure(
          "FramePlan",
          s"frame index $frameIndex out of range for signal length ${signal.length}"
        )
      )
    else
      val start = frameStart(frameIndex)
      val out = DVecBuilder.zeros(frameLength)
      var i = 0
      while i < frameLength do
        out(i) = SampleAccess(signal, start + i, boundary) * window.taps(i)
        i += 1
      Right(out.result())

  /** Overlap-add `frame` (already synthesis-windowed) into `dest` at frame index. */
  def overlapAdd(dest: Array[Double], frame: DVec, frameIndex: Int): Unit =
    val start = frameStart(frameIndex)
    var i = 0
    while i < frame.length do
      val idx = start + i
      if idx >= 0 && idx < dest.length then dest(idx) += frame(i)
      i += 1

object FramePlan:
  def apply(
      frameLength: Int,
      hop: Int,
      alignment: FrameAlignment,
      boundary: Boundary = Boundary.Zero
  ): Either[SignalError, FramePlan] =
    if frameLength <= 0 then Left(SignalError.InvalidInputLength(frameLength))
    else if hop <= 0 then Left(SignalError.InvalidInputLength(hop))
    else Right(new FramePlan(frameLength, hop, alignment, boundary))
