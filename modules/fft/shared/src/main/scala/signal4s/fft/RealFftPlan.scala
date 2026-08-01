package signal4s.fft

import gale.linalg.{DVec, DVecBuilder}
import signal4s.{SampleRate, SignalError}
import signal4s.fft.internal.FftEngine

/** Immutable real-input FFT plan. Forward yields [[RealSpectrum]] with
  * `n/2+1` complex bins (NumPy `rfft` layout).
  */
final class RealFftPlan private (
    val length: Int,
    val normalization: FftNormalization,
    val sampleRate: Option[SampleRate]
):
  def binCount: Int = length / 2 + 1

  def newWorkspace(): FftWorkspace =
    val packLength = if usesPackedReal then length / 2 else length
    FftWorkspace(
      Array.ofDim[Double](length),
      Array.ofDim[Double](length),
      Array.ofDim[Double](packLength),
      Array.ofDim[Double](packLength)
    )

  def forward(x: DVec): Either[SignalError, RealSpectrum] =
    val ws = newWorkspace()
    forwardInto(x, ws).flatMap(_ => spectrumFrom(ws))

  def inverse(spectrum: RealSpectrum): Either[SignalError, DVec] =
    val ws = newWorkspace()
    inverseInto(spectrum, ws).map(_ => realOut(ws))

  def forwardInto(input: DVec, workspace: FftWorkspace): Either[SignalError, Unit] =
    if input.length != length then Left(SignalError.LengthMismatch(length, input.length))
    else if workspace.length != length then
      Left(SignalError.LengthMismatch(length, workspace.length))
    else
      var i = 0
      while i < length do
        workspace.re(i) = input(i)
        workspace.im(i) = 0.0
        i += 1
      if usesPackedReal then
        FftEngine.forwardRealOnesided(
          workspace.re,
          workspace.im,
          workspace.scratchRe,
          workspace.scratchIm
        )
      else FftEngine.forward(workspace.re, workspace.im)
      applyScale(workspace, forward = true)
      Right(())

  def inverseInto(
      spectrum: RealSpectrum,
      workspace: FftWorkspace
  ): Either[SignalError, Unit] =
    if spectrum.sourceLength != length then
      Left(SignalError.LengthMismatch(length, spectrum.sourceLength))
    else if spectrum.bins.length != binCount then
      Left(SignalError.LengthMismatch(binCount, spectrum.bins.length))
    else if workspace.length != length then
      Left(SignalError.LengthMismatch(length, workspace.length))
    else
      workspace.re(0) = spectrum.bins.real(0)
      workspace.im(0) = spectrum.bins.imaginary(0)
      var k = 1
      while k < binCount do
        val re = spectrum.bins.real(k)
        val im = spectrum.bins.imaginary(k)
        workspace.re(k) = re
        workspace.im(k) = im
        if !usesPackedReal && k < length - k then
          workspace.re(length - k) = re
          workspace.im(length - k) = -im
        k += 1
      if usesPackedReal then
        FftEngine.inverseReal(
          workspace.re,
          workspace.im,
          workspace.scratchRe,
          workspace.scratchIm
        )
      else
        // Hermitian expansion above supplies the full spectrum for arbitrary lengths.
        FftEngine.inverse(workspace.re, workspace.im)
      applyScale(workspace, forward = false)
      Right(())

  def spectrumFrom(workspace: FftWorkspace): Either[SignalError, RealSpectrum] =
    if workspace.length != length then
      Left(SignalError.LengthMismatch(length, workspace.length))
    else
      val re = DVecBuilder.zeros(binCount)
      val im = DVecBuilder.zeros(binCount)
      var k = 0
      while k < binCount do
        re(k) = workspace.re(k)
        im(k) = workspace.im(k)
        k += 1
      val bins = ComplexVector.unsafe(re.result(), im.result())
      val axis =
        sampleRate match
          case Some(fs) => FrequencyAxis.realFft(length, fs)
          case None     => FrequencyAxis.normalizedRealFft(length)
      axis.flatMap(a => RealSpectrum(bins, a, length, normalization))

  private def realOut(ws: FftWorkspace): DVec =
    val out = DVecBuilder.zeros(length)
    var i = 0
    while i < length do
      out(i) = ws.re(i)
      i += 1
    out.result()

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

  private def usesPackedReal: Boolean =
    FftEngine.supportsRealLength(length)

object RealFftPlan:
  def apply(
      length: Int,
      normalization: FftNormalization
  ): Either[SignalError, RealFftPlan] =
    if length <= 0 then Left(SignalError.InvalidInputLength(length))
    else Right(new RealFftPlan(length, normalization, None))

  def apply(
      length: Int,
      normalization: FftNormalization,
      sampleRate: SampleRate
  ): Either[SignalError, RealFftPlan] =
    if length <= 0 then Left(SignalError.InvalidInputLength(length))
    else Right(new RealFftPlan(length, normalization, Some(sampleRate)))
