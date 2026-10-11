package signal4s.fft

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*

/** Numerical array payloads only; excludes objects, references, allocator, GC and RSS. */
final case class PeriodogramResources(windowBytes: BigInt, workspaceBytes: BigInt,
    frameScratchBytes: BigInt, ownedResultBytes: BigInt, fftWorkingLength: Int):
  def peakAdditionalBytes: BigInt = workspaceBytes + frameScratchBytes + ownedResultBytes

/** Single-owner reusable arrays/tables. Shape/owner refusals occur before modification. */
final class PeriodogramWorkspace private[fft] (private[fft] val owner: BoundedPeriodogramPlan):
  private[fft] val fft = owner.fftLayout.newWorkspace()
  private[fft] val re = fft.re
  private[fft] val im = fft.im
  private[fft] val tables = fft.tables
  private[fft] val chirpRe = fft.chirpRe
  private[fft] val chirpIm = fft.chirpIm
  private[fft] val kernelRe = fft.kernelRe
  private[fft] val kernelIm = fft.kernelIm

/** One full frame, explicit portable FFT ownership, no native/global numerical caches.
  * Radix-2 is reused from the engine; arbitrary lengths use owned Bluestein arrays.
  * Resource inspection and construction allocate no FFT workspace or tables.
  */
final class BoundedPeriodogramPlan private (val window: Window, val sampleRate: SampleRate,
    val nfft: Int, val detrend: Detrend, val scaling: SpectralScaling,
    val resources: PeriodogramResources, private val divisor: Double,
    private val windowPower: Double, private val coherentPower: Double):
  private[fft] val fftLayout = OwnedFrameFft(nfft).orThrow
  def binCount: Int = nfft / 2 + 1
  def newWorkspace(): PeriodogramWorkspace = new PeriodogramWorkspace(this)
  def estimate(input: DVec): Either[SignalError, WelchResult] = estimateInto(input, newWorkspace())

  def estimateInto(input: DVec, workspace: PeriodogramWorkspace): Either[SignalError, WelchResult] =
    if workspace.owner ne this then Left(SignalError.NumericalFailure("periodogram.workspace", "workspace belongs to a different plan"))
    else SpectralFrame.prepare(input, window, detrend).flatMap { case (frame, scale) =>
      workspace.fft.forward(frame)
      val out = DVecBuilder.zeros(binCount)
      var i = 0
      var error: Option[SignalError] = None
      while i < binCount && error.isEmpty do
        val re = workspace.fft.real(i)
        val im = workspace.fft.imaginary(i)
        val factor = if i == 0 || (nfft % 2 == 0 && i == binCount - 1) then 1.0 else 2.0
        SpectralFrame.power(re, im, scale, divisor, factor) match
          case Left(e) => error = Some(e)
          case Right(value) => out(i) = value
        i += 1
      error match
        case Some(e) => Left(e)
        case None => FrequencyAxis.realFft(nfft, sampleRate).map(axis => WelchResult(axis, out.result(),
          scaling, SpectralSides.Onesided, 1, windowPower, coherentPower, detrend, AverageMethod.Mean,
          None, List("portable owned FFT", "one full frame", "dof=unspecified (not fabricated)")))
    }

object BoundedPeriodogramPlan:
  /** Bounds the Bluestein working length and all Int-indexed arrays before allocation. */
  val MaxFftLength: Int = OwnedFrameFft.MaxFftLength
  def apply(window: Window, sampleRate: SampleRate, nfft: Int,
      detrend: Detrend = Detrend.Mean, scaling: SpectralScaling = SpectralScaling.Density): Either[SignalError, BoundedPeriodogramPlan] =
    if nfft < window.length || nfft > MaxFftLength then
      Left(SignalError.NumericalFailure("periodogram.shape", "FFT length outside window/bounded capacity"))
    else
      FrequencyAxis.realFft(nfft, sampleRate).flatMap { _ =>
        WelchPlan(window, window.length, sampleRate, nfft, detrend = detrend, scaling = scaling)
      }.flatMap { _ =>
        var sumW = 0.0
        var sumW2 = 0.0
        var i = 0
        while i < window.length do
          sumW += window.taps(i)
          sumW2 += window.taps(i) * window.taps(i)
          i += 1
        val layout = OwnedFrameFft(nfft).orThrow
        val r = PeriodogramResources(BigInt(window.length) * 8, layout.workspaceBytes,
          BigInt(window.length) * 64 + gale.numeric.ExactSum.resources.stateBytes * 2 +
            gale.numeric.ExactSum.resources.ratioScratchBytes,
          BigInt(nfft / 2 + 1) * 16, layout.workingLength)
        val divisor = scaling match
          case SpectralScaling.Density => sampleRate.hertz * sumW2
          case SpectralScaling.Spectrum => sumW * sumW
        Right(new BoundedPeriodogramPlan(window, sampleRate, nfft, detrend, scaling, r, divisor, sumW2, sumW * sumW))
      }
