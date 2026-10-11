package signal4s.fft

import gale.linalg.{DVec, DVecBuilder}
import gale.numeric.ExactSum
import signal4s.*

/** Complex positive-frequency coefficients, without interior-bin folding. */
enum StftScaling:
  case Raw, Magnitude, Density

/** Primitive/limb payload capacities; excludes objects, GC and RSS. */
final case class StftResources(windowBytes: BigInt, workspaceBytes: BigInt,
    frameScratchBytes: BigInt, ownedResultBytes: BigInt,
    constructionAdditionalBytes: BigInt, fftWorkingLength: Int):
  def peakBytes: BigInt = windowBytes + workspaceBytes + frameScratchBytes + ownedResultBytes

/** Detached full-frame analysis result; its phase origin is the input's first sample. */
final class StftFrame private[fft] (val bins: ComplexVector, val frequencies: FrequencyAxis,
    val nfft: Int, val frameLength: Int, val scaling: StftScaling, val detrend: Detrend):
  val sides: SpectralSides = SpectralSides.Onesided
  val phaseOrigin: String = "frame_start; raw rFFT; phase_shift=None"

/** Mutable single-owner scratch. A frame result never borrows these arrays. */
final class StftWorkspace private[fft] (private[fft] val owner: BoundedStftPlan):
  private[fft] val fft = owner.fftLayout.newWorkspace()

/** Analysis only: one complete observation frame, one owned portable FFT workspace.
  * There is no hop, padding boundary, frame collection or synthesis requirement.
  */
final class BoundedStftPlan private (val window: Window, val sampleRate: SampleRate,
    val nfft: Int, val detrend: Detrend, val scaling: StftScaling,
    val resources: StftResources, private val normalizedWindow: Window,
    private val windowScale: Double, private val coherent: Double,
    private val energy: Double, private[fft] val fftLayout: OwnedFrameFft):
  val frameLength: Int = window.length
  val binCount: Int = nfft / 2 + 1

  /** Budget includes retained plan taps and one frame/result, before FFT allocation. */
  def newWorkspace(numericBudgetBytes: BigInt): Either[SignalError, StftWorkspace] =
    if numericBudgetBytes < resources.peakBytes then
      Left(SignalError.NumericalFailure("stft.budget", "numerical capacity below declared plan/frame/workspace/result peak"))
    else Right(new StftWorkspace(this))

  def analyze(input: DVec, numericBudgetBytes: BigInt): Either[SignalError, StftFrame] =
    if input.length != frameLength then Left(SignalError.LengthMismatch(frameLength, input.length))
    else if !input.toSeq.forall(_.isFinite) then
      Left(SignalError.NumericalFailure("stft.frame", "nonfinite input"))
    else newWorkspace(numericBudgetBytes).flatMap(analyzeInto(input, _))

  def analyzeInto(input: DVec, workspace: StftWorkspace): Either[SignalError, StftFrame] =
    if workspace.owner ne this then
      Left(SignalError.NumericalFailure("stft.workspace", "workspace belongs to a different plan"))
    else SpectralFrame.prepare(input, normalizedWindow, detrend).flatMap { case (frame, scale) =>
      workspace.fft.forward(frame)
      val real = DVecBuilder.zeros(binCount)
      val imaginary = DVecBuilder.zeros(binCount)
      var k = 0
      var error: Option[SignalError] = None
      while k < binCount && error.isEmpty do
        val re = restore(workspace.fft.real(k), scale)
        val im = restore(workspace.fft.imaginary(k), scale)
        (re, im) match
          case (Right(r), Right(i)) => real(k) = r; imaginary(k) = i
          case (Left(e), _) => error = Some(e)
          case (_, Left(e)) => error = Some(e)
        k += 1
      error match
        case Some(e) => Left(e)
        case None => FrequencyAxis.realFft(nfft, sampleRate).map(axis =>
          new StftFrame(ComplexVector.unsafe(real.result(), imaginary.result()), axis,
            nfft, frameLength, scaling, detrend))
    }

  private def restore(component: Double, scale: Double): Either[SignalError, Double] =
    if !component.isFinite then Left(SignalError.NumericalFailure("stft.fft", "nonfinite FFT coefficient"))
    else if component == 0 || scale == 0 || windowScale == 0 then Right(component * 0.0)
    else
      // Restore only after normalization; no overflowing raw amplitude/divisor products.
      var mantissa = 1.0
      var exponent = 0
      def combine(value: Double, divide: Boolean): Unit =
        val initial = java.lang.Math.getExponent(math.abs(value))
        val e = if initial == -1023 then java.lang.Math.getExponent(java.lang.Math.scalb(math.abs(value), 54)) - 54 else initial
        val m = java.lang.Math.scalb(value, -e)
        mantissa = if divide then mantissa / m else mantissa * m
        exponent += (if divide then -e else e)
      combine(component, false)
      combine(scale, false)
      scaling match
        case StftScaling.Raw => combine(windowScale, false)
        case StftScaling.Magnitude => combine(coherent, true)
        case StftScaling.Density =>
          combine(math.sqrt(sampleRate.hertz), true)
          combine(math.sqrt(energy), true)
      val value = java.lang.Math.scalb(mantissa, exponent)
      if value.isFinite then Right(value)
      else Left(SignalError.NumericalFailure("stft.coefficient", "normalized coefficient exceeds finite capacity"))

object BoundedStftPlan:
  val MaxFftLength: Int = OwnedFrameFft.MaxFftLength

  /** Shape-only admission. Does not allocate window/frame/FFT buffers. */
  def inspect(frameLength: Int, nfft: Int): Either[SignalError, StftResources] =
    if frameLength <= 0 || nfft < frameLength then
      Left(SignalError.NumericalFailure("stft.shape", "positive complete frame and FFT at least its length required"))
    else OwnedFrameFft(nfft).map { layout =>
      val exact = ExactSum.resources
      val scratch = BigInt(frameLength) * 64 + BigInt(exact.stateBytes) * 2 + exact.ratioScratchBytes
      StftResources(BigInt(frameLength) * 16, layout.workspaceBytes, scratch,
        BigInt(nfft / 2 + 1) * 32, scratch + BigInt(frameLength) * 16, layout.workingLength)
    }

  /** Admit construction and the one-frame run before copying any window taps. */
  def withinBudget(window: Window, sampleRate: SampleRate, nfft: Int,
      numericBudgetBytes: BigInt, detrend: Detrend = Detrend.None,
      scaling: StftScaling = StftScaling.Raw): Either[SignalError, BoundedStftPlan] =
    inspect(window.length, nfft).flatMap { r =>
      if numericBudgetBytes < r.peakBytes.max(r.constructionAdditionalBytes + r.windowBytes) then
        Left(SignalError.NumericalFailure("stft.budget", "capacity below construction or one-frame peak"))
      else apply(window, sampleRate, nfft, detrend, scaling)
    }

  def apply(window: Window, sampleRate: SampleRate, nfft: Int,
      detrend: Detrend = Detrend.None, scaling: StftScaling = StftScaling.Raw): Either[SignalError, BoundedStftPlan] =
    for
      resources <- inspect(window.length, nfft)
      _ <- if sampleRate.hertz.isFinite && sampleRate.hertz > 0 && sampleRate.hertz / nfft > 0 then Right(())
        else Left(SignalError.NumericalFailure("stft.frequency", "positive finite rate and representable positive FFT spacing required"))
      owned <- Window.fromTaps(window.taps.copy, window.convention)
      maximum = owned.taps.toSeq.map(math.abs).max
      normalized <- Window.fromTaps(DVec.tabulate(owned.length)(i =>
        if maximum == 0 then 0.0 else owned.taps(i) / maximum), owned.convention)
      coherent <- sum(normalized.taps, square = false)
      energy <- sum(normalized.taps, square = true)
      _ <- scaling match
        case StftScaling.Magnitude if coherent == 0 =>
          Left(SignalError.NumericalFailure("stft.normalization", "magnitude normalization requires nonzero coherent window sum"))
        case StftScaling.Density if energy <= 0 =>
          Left(SignalError.NumericalFailure("stft.normalization", "density normalization requires positive window energy"))
        case _ => Right(())
      layout <- OwnedFrameFft(nfft)
    yield new BoundedStftPlan(owned, sampleRate, nfft, detrend, scaling, resources,
      normalized, maximum, math.abs(coherent), energy, layout)

  private def sum(values: DVec, square: Boolean): Either[SignalError, Double] =
    val accumulator = ExactSum.zero()
    var i = 0
    while i < values.length do
      val v = values(i)
      accumulator.add(if square then v * v else v) match
        case Left(e) => return Left(SignalError.NumericalFailure("stft.window", e.message))
        case Right(_) => ()
      i += 1
    Right(accumulator.value)
