package signal4s.fft

import gale.linalg.DVecBuilder
import signal4s.SignalError
import signal4s.fft.internal.FftEngine

/** Immutable, thread-safe complex FFT plan of fixed length. */
final class FftPlan private (
    val length: Int,
    val normalization: FftNormalization
):
  def newWorkspace(): FftWorkspace =
    FftWorkspace(
      Array.ofDim[Double](length),
      Array.ofDim[Double](length),
      Array.ofDim[Double](length),
      Array.ofDim[Double](length)
    )

  def forward(x: ComplexVector): Either[SignalError, ComplexVector] =
    val ws = newWorkspace()
    forwardInto(x, ws).map(_ => copyOut(ws))

  def inverse(x: ComplexVector): Either[SignalError, ComplexVector] =
    val ws = newWorkspace()
    inverseInto(x, ws).map(_ => copyOut(ws))

  def forwardInto(
      input: ComplexVector,
      workspace: FftWorkspace
  ): Either[SignalError, Unit] =
    if input.length != length then Left(SignalError.LengthMismatch(length, input.length))
    else if workspace.length != length then
      Left(SignalError.LengthMismatch(length, workspace.length))
    else
      copyIn(input, workspace)
      FftEngine.forward(workspace.re, workspace.im)
      applyScale(workspace, forward = true)
      Right(())

  def inverseInto(
      input: ComplexVector,
      workspace: FftWorkspace
  ): Either[SignalError, Unit] =
    if input.length != length then Left(SignalError.LengthMismatch(length, input.length))
    else if workspace.length != length then
      Left(SignalError.LengthMismatch(length, workspace.length))
    else
      copyIn(input, workspace)
      FftEngine.inverse(workspace.re, workspace.im)
      applyScale(workspace, forward = false)
      Right(())

  /** Read workspace buffers into a new [[ComplexVector]]. */
  def resultFrom(workspace: FftWorkspace): Either[SignalError, ComplexVector] =
    if workspace.length != length then
      Left(SignalError.LengthMismatch(length, workspace.length))
    else Right(copyOut(workspace))

  private def copyIn(input: ComplexVector, ws: FftWorkspace): Unit =
    var i = 0
    while i < length do
      ws.re(i) = input.real(i)
      ws.im(i) = input.imaginary(i)
      i += 1

  private def copyOut(ws: FftWorkspace): ComplexVector =
    val re = DVecBuilder.zeros(length)
    val im = DVecBuilder.zeros(length)
    var i = 0
    while i < length do
      re(i) = ws.re(i)
      im(i) = ws.im(i)
      i += 1
    ComplexVector.unsafe(re.result(), im.result())

  private def applyScale(ws: FftWorkspace, forward: Boolean): Unit =
    val scale =
      normalization match
        case FftNormalization.Backward =>
          if forward then 1.0 else 1.0 / length.toDouble
        case FftNormalization.Forward =>
          if forward then 1.0 / length.toDouble else 1.0
        case FftNormalization.Orthonormal =>
          1.0 / math.sqrt(length.toDouble)
    if scale != 1.0 then
      var i = 0
      while i < length do
        ws.re(i) *= scale
        ws.im(i) *= scale
        i += 1

object FftPlan:
  // Touching the FFT plan API installs convolution hooks into core.
  FftBackend.ensureInstalled()

  def apply(
      length: Int,
      normalization: FftNormalization = FftNormalization.Backward
  ): Either[SignalError, FftPlan] =
    if length <= 0 then Left(SignalError.InvalidInputLength(length))
    else Right(new FftPlan(length, normalization))
