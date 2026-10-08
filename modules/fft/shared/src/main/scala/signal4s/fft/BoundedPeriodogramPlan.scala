package signal4s.fft

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*
import signal4s.fft.internal.{FftEngine, Radix2Tables}

/** Numerical array payloads only; excludes objects, references, allocator, GC and RSS. */
final case class PeriodogramResources(windowBytes: BigInt, workspaceBytes: BigInt,
    frameScratchBytes: BigInt, ownedResultBytes: BigInt, fftWorkingLength: Int):
  def peakAdditionalBytes: BigInt = workspaceBytes + frameScratchBytes + ownedResultBytes

/** Single-owner reusable arrays/tables. Shape/owner refusals occur before modification. */
final class PeriodogramWorkspace private[fft] (private[fft] val owner: BoundedPeriodogramPlan):
  private[fft] val re = Array.ofDim[Double](owner.resources.fftWorkingLength)
  private[fft] val im = Array.ofDim[Double](owner.resources.fftWorkingLength)
  private[fft] val tables = Radix2Tables.owned(owner.resources.fftWorkingLength)
  private[fft] val chirpRe = if owner.isRadix2 then Array.emptyDoubleArray else Array.ofDim[Double](owner.nfft)
  private[fft] val chirpIm = if owner.isRadix2 then Array.emptyDoubleArray else Array.ofDim[Double](owner.nfft)
  private[fft] val kernelRe = if owner.isRadix2 then Array.emptyDoubleArray else Array.ofDim[Double](re.length)
  private[fft] val kernelIm = if owner.isRadix2 then Array.emptyDoubleArray else Array.ofDim[Double](re.length)
  if !owner.isRadix2 then
    var i = 0
    while i < owner.nfft do
      val angle = -math.Pi * ((i.toLong * i) % (2L * owner.nfft)).toDouble / owner.nfft.toDouble
      chirpRe(i) = math.cos(angle)
      chirpIm(i) = math.sin(angle)
      kernelRe(i) = chirpRe(i)
      kernelIm(i) = -chirpIm(i)
      if i > 0 then
        kernelRe(re.length - i) = chirpRe(i)
        kernelIm(re.length - i) = -chirpIm(i)
      i += 1
    FftEngine.radix2Owned(kernelRe, kernelIm, false, tables)

/** One full frame, explicit portable FFT ownership, no native/global numerical caches.
  * Radix-2 is reused from the engine; arbitrary lengths use owned Bluestein arrays.
  * Resource inspection and construction allocate no FFT workspace or tables.
  */
final class BoundedPeriodogramPlan private (val window: Window, val sampleRate: SampleRate,
    val nfft: Int, val detrend: Detrend, val scaling: SpectralScaling,
    val resources: PeriodogramResources, private val divisor: Double,
    private val windowPower: Double, private val coherentPower: Double):
  private[fft] val isRadix2 = (nfft & (nfft - 1)) == 0
  def binCount: Int = nfft / 2 + 1
  def newWorkspace(): PeriodogramWorkspace = new PeriodogramWorkspace(this)
  def estimate(input: DVec): Either[SignalError, WelchResult] = estimateInto(input, newWorkspace())

  def estimateInto(input: DVec, workspace: PeriodogramWorkspace): Either[SignalError, WelchResult] =
    if workspace.owner ne this then Left(SignalError.NumericalFailure("periodogram.workspace", "workspace belongs to a different plan"))
    else SpectralFrame.prepare(input, window, detrend).flatMap { case (frame, scale) =>
      java.util.Arrays.fill(workspace.re, 0.0)
      java.util.Arrays.fill(workspace.im, 0.0)
      var i = 0
      while i < frame.length do
        if isRadix2 then workspace.re(i) = frame(i)
        else
          workspace.re(i) = frame(i) * workspace.chirpRe(i)
          workspace.im(i) = frame(i) * workspace.chirpIm(i)
        i += 1
      FftEngine.radix2Owned(workspace.re, workspace.im, false, workspace.tables)
      if !isRadix2 then
        i = 0
        while i < workspace.re.length do
          val ar = workspace.re(i)
          val ai = workspace.im(i)
          workspace.re(i) = ar * workspace.kernelRe(i) - ai * workspace.kernelIm(i)
          workspace.im(i) = ar * workspace.kernelIm(i) + ai * workspace.kernelRe(i)
          i += 1
        FftEngine.radix2Owned(workspace.re, workspace.im, true, workspace.tables)
      val out = DVecBuilder.zeros(binCount)
      i = 0
      var error: Option[SignalError] = None
      while i < binCount && error.isEmpty do
        val re = if isRadix2 then workspace.re(i) else
          (workspace.re(i) * workspace.chirpRe(i) - workspace.im(i) * workspace.chirpIm(i)) / workspace.re.length
        val im = if isRadix2 then workspace.im(i) else
          (workspace.re(i) * workspace.chirpIm(i) + workspace.im(i) * workspace.chirpRe(i)) / workspace.re.length
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
  val MaxFftLength: Int = 1 << 28
  def apply(window: Window, sampleRate: SampleRate, nfft: Int,
      detrend: Detrend = Detrend.Mean, scaling: SpectralScaling = SpectralScaling.Density): Either[SignalError, BoundedPeriodogramPlan] =
    if nfft < window.length || nfft > MaxFftLength then
      Left(SignalError.NumericalFailure("periodogram.shape", "FFT length outside window/bounded capacity"))
    else
      WelchPlan(window, window.length, sampleRate, nfft, detrend = detrend, scaling = scaling).flatMap { _ =>
        var sumW = 0.0
        var sumW2 = 0.0
        var i = 0
        while i < window.length do
          sumW += window.taps(i)
          sumW2 += window.taps(i) * window.taps(i)
          i += 1
        val radix2 = (nfft & (nfft - 1)) == 0
        var working = 1
        val minimum = if radix2 then nfft.toLong else 2L * nfft - 1
        while working.toLong < minimum do working <<= 1
        val n = BigInt(working)
        // Both forward/inverse tables: bitrev Int[n] + four Double stages totaling n-1 taps each.
        val tables = 4 * n + 32 * (n - 1)
        val arrays = 16 * n + (if radix2 then BigInt(0) else 16 * n + 16 * BigInt(nfft))
        val r = PeriodogramResources(BigInt(window.length) * 8, tables + arrays,
          BigInt(window.length) * 64 + gale.numeric.ExactSum.resources.stateBytes * 2 +
            gale.numeric.ExactSum.resources.ratioScratchBytes,
          BigInt(nfft / 2 + 1) * 16, working)
        val divisor = scaling match
          case SpectralScaling.Density => sampleRate.hertz * sumW2
          case SpectralScaling.Spectrum => sumW * sumW
        Right(new BoundedPeriodogramPlan(window, sampleRate, nfft, detrend, scaling, r, divisor, sumW2, sumW * sumW))
      }
