package signal4s.fft

import gale.linalg.{DVec, DVecBuilder}
import gale.numeric.ExactSum
import signal4s.{SampleRate, SignalError}

/** Primitive payload bound, excluding caller inputs, object headers, GC and RSS. */
final case class CsdResources(stateBytes: BigInt, scratchBytes: BigInt,
    outputBytes: BigInt, constructionBytes: BigInt):
  def peakBytes: BigInt = (stateBytes + scratchBytes + outputBytes).max(constructionBytes)

/** Owned full matrix in ((row * channels + column) * bins + frequency) order.
  * S(i,j) = mean(conj(X_i) * X_j), one-sided density. No degrees-of-freedom claim.
  */
sealed trait CsdEstimate:
  def channels: Int
  def binCount: Int
  def real: DVec
  def imaginary: DVec
  def acceptedWindows: Long
  def sourceFrequencies: FrequencyAxis
  def frequencyGroups: Vector[Vector[Int]]
  final def apply(row: Int, column: Int, bin: Int): Complex =
    require(row >= 0 && row < channels && column >= 0 && column < channels && bin >= 0 && bin < binCount)
    val index = (row * channels + column) * binCount + bin
    Complex(real(index), imaginary(index))

final class CsdMatrix private[fft] (val channels: Int, val binCount: Int,
    val real: DVec, val imaginary: DVec, val acceptedWindows: Long,
    val frequencies: FrequencyAxis) extends CsdEstimate:
  def sourceFrequencies: FrequencyAxis = frequencies
  val frequencyGroups: Vector[Vector[Int]] = Vector.tabulate(binCount)(i => Vector(i))
  /** Common nonnegative averaging preserves density units, rather than integrating power. */
  def averageBins(groups: Vector[Vector[Int]], numericBudgetBytes: BigInt): Either[SignalError, CsdGroupedMatrix] =
    for
      bytes <- CsdGroupedMatrix.inspect(channels, binCount, groups)
      _ <- if bytes <= numericBudgetBytes then Right(()) else Left(BoundedCsdPlan.failure("grouped CSD output exceeds budget"))
    yield
      val size = channels * channels * groups.size
      val outputReal = DVecBuilder.zeros(size)
      val outputImaginary = DVecBuilder.zeros(size)
      for row <- 0 until channels; column <- row until channels; group <- groups.indices do
        val re = ExactSum.zero(); val im = ExactSum.zero(); val count = ExactSum.zero()
        groups(group).foreach: bin =>
          re.add(apply(row, column, bin).real).fold(e => throw new IllegalStateException(e.toString), identity)
          im.add(apply(row, column, bin).imaginary).fold(e => throw new IllegalStateException(e.toString), identity)
        count.add(groups(group).size.toDouble)
        val realMean = re.ratio(count).fold(e => throw new IllegalStateException(e.toString), identity)
        val imaginaryMean = if row == column then 0.0 else im.ratio(count).fold(e => throw new IllegalStateException(e.toString), identity)
        val upper = (row * channels + column) * groups.size + group
        val lower = (column * channels + row) * groups.size + group
        outputReal(upper) = realMean; outputReal(lower) = realMean
        outputImaginary(upper) = imaginaryMean
        outputImaginary(lower) = if row == column then 0.0 else -imaginaryMean
      new CsdGroupedMatrix(channels, outputReal.result(), outputImaginary.result(), acceptedWindows, frequencies, groups)

/** Owned arithmetic group means; source bins and their original axis remain explicit. */
final class CsdGroupedMatrix private[fft] (val channels: Int, val real: DVec, val imaginary: DVec,
    val acceptedWindows: Long, val sourceFrequencies: FrequencyAxis,
    val frequencyGroups: Vector[Vector[Int]]) extends CsdEstimate:
  val binCount: Int = frequencyGroups.size

object CsdGroupedMatrix:
  /** Includes the retained full source matrix and builder/result overlap. */
  def inspect(channels: Int, bins: Int, groups: Vector[Vector[Int]]): Either[SignalError, BigInt] =
    var previous = -1
    val ordered = groups.nonEmpty && groups.forall: group =>
      group.nonEmpty && group.forall: bin =>
        val valid = bin > previous && bin < bins
        previous = bin
        valid
    val full = BigInt(channels) * channels * bins
    val result = BigInt(channels) * channels * groups.size
    if channels <= 0 || bins <= 0 || !ordered || full > Int.MaxValue || result > Int.MaxValue then
      Left(BoundedCsdPlan.failure("nonempty ordered disjoint frequency groups inside bounded CSD shape required"))
    else Right(full * 16 + BigInt(bins) * 4 + result * 32 + BigInt(groups.iterator.map(_.size).sum) * 4 +
      ExactSum.resources.stateBytes * 3 + ExactSum.resources.ratioScratchBytes)

/** Common-window estimator; the caller determines which complete windows exist.
  * A single invalid row refuses the entire frame before any spectral sum changes.
  */
final class BoundedCsdPlan private (val channels: Int, val stft: BoundedStftPlan,
    val resources: CsdResources):
  val binCount: Int = stft.binCount
  private[fft] val cells: Int = (channels.toLong * (channels + 1) / 2 * binCount).toInt
  def newAccumulator(numericBudgetBytes: BigInt): Either[SignalError, CsdAccumulator] =
    if numericBudgetBytes < resources.peakBytes then Left(BoundedCsdPlan.failure("budget below declared CSD peak"))
    else stft.newWorkspace(numericBudgetBytes).map(w => new CsdAccumulator(this, w))

object BoundedCsdPlan:
  private[fft] def failure(message: String): SignalError = SignalError.NumericalFailure("csd", message)
  val MaxWindows: Long = math.min(ExactSum.MaxTerms, 9007199254740991L)

  def inspect(channels: Int, frameLength: Int, nfft: Int): Either[SignalError, CsdResources] =
    val bins = BigInt(nfft) / 2 + 1
    val full = BigInt(channels) * channels * bins
    val triangle = BigInt(channels) * (channels + 1) / 2 * bins
    if channels <= 0 || full > Int.MaxValue || triangle > Int.MaxValue then Left(failure("positive bounded matrix shape required"))
    else BoundedStftPlan.inspect(frameLength, nfft).map: fft =>
      val exact = ExactSum.resources
      val state = triangle * 2 * exact.stateBytes
      val output = full * 16 * 2 + bins * 4 // Builders, detached values and original bin membership.
      val scratch = triangle * 16 + BigInt(channels) * bins * 32 + fft.peakBytes +
        exact.ratioScratchBytes + exact.stateBytes
      CsdResources(state, scratch, output, fft.constructionAdditionalBytes + fft.windowBytes)

  def withinBudget(channels: Int, window: Window, sampleRate: SampleRate, nfft: Int,
      numericBudgetBytes: BigInt, detrend: Detrend = Detrend.None): Either[SignalError, BoundedCsdPlan] =
    for
      resource <- inspect(channels, window.length, nfft)
      _ <- if resource.peakBytes <= numericBudgetBytes then Right(()) else Left(failure("budget below declared CSD peak"))
      stft <- BoundedStftPlan.withinBudget(window, sampleRate, nfft, numericBudgetBytes, detrend, StftScaling.Density)
    yield new BoundedCsdPlan(channels, stft, resource)

/** Ordinary and multitaper estimators expose the same common-frame transaction. */
trait CsdFrameAccumulator:
  def acceptedWindows: Long
  def accumulateFrame(rows: IndexedSeq[DVec]): Either[SignalError,Unit]
  def finishOwned(): Either[SignalError,Option[CsdMatrix]]

/** Single owner; frequency states are fixed-size and a frame commits atomically. */
final class CsdAccumulator private[fft] (val plan: BoundedCsdPlan, workspace: StftWorkspace) extends CsdFrameAccumulator:
  private val real = Array.fill(plan.cells)(ExactSum.zero())
  private val imaginary = Array.fill(plan.cells)(ExactSum.zero())
  private val stagedReal = new Array[Double](plan.cells)
  private val stagedImaginary = new Array[Double](plan.cells)
  private var count = 0L
  def acceptedWindows: Long = count

  def accumulateFrame(rows: IndexedSeq[DVec]): Either[SignalError, Unit] =
    if rows.length != plan.channels then return Left(BoundedCsdPlan.failure("channel count differs from plan"))
    if count == BoundedCsdPlan.MaxWindows then return Left(BoundedCsdPlan.failure("exact window count capacity exceeded"))
    val spectra = new Array[ComplexVector](plan.channels)
    var row = 0
    while row < plan.channels do
      plan.stft.analyzeInto(rows(row), workspace) match
        case Left(error) => return Left(error)
        case Right(frame) => spectra(row) = frame.bins
      row += 1
    CrossPeriodogram.write(spectra,plan.stft.nfft,stagedReal,stagedImaginary).flatMap(_=>accumulateSpectrum(stagedReal,stagedImaginary))

  /** Internal upper-triangle frame contribution, fully preflighted before commit. */
  private[fft] def accumulateSpectrum(frameReal: Array[Double],frameImaginary: Array[Double]): Either[SignalError,Unit] =
    if count==BoundedCsdPlan.MaxWindows then return Left(BoundedCsdPlan.failure("exact window count capacity exceeded"))
    if frameReal.length!=plan.cells || frameImaginary.length!=plan.cells ||
      !frameReal.forall(_.isFinite) || !frameImaginary.forall(_.isFinite) then
      return Left(BoundedCsdPlan.failure("invalid common-frame spectrum shape or finite capacity"))
    // Every state receives at most one term per frame. The global exact count
    // preflight proves all add capacities before this first observable update.
    var index = 0
    while index < plan.cells do
      real(index).add(frameReal(index)).fold(e => throw new IllegalStateException(e.toString), identity)
      imaginary(index).add(frameImaginary(index)).fold(e => throw new IllegalStateException(e.toString), identity)
      index += 1
    count += 1
    Right(())

  def merge(other: CsdAccumulator): Either[SignalError, Unit] =
    if this eq other then Left(BoundedCsdPlan.failure("self merge duplicates observations"))
    else if plan ne other.plan then Left(BoundedCsdPlan.failure("different plan identity"))
    else if other.count > BoundedCsdPlan.MaxWindows - count then Left(BoundedCsdPlan.failure("exact window count capacity exceeded"))
    else
      var index = 0
      while index < plan.cells do
        real(index).addAll(other.real(index)).fold(e => throw new IllegalStateException(e.toString), identity)
        imaginary(index).addAll(other.imaginary(index)).fold(e => throw new IllegalStateException(e.toString), identity)
        index += 1
      count += other.count
      Right(())

  def finishOwned(): Either[SignalError, Option[CsdMatrix]] =
    if count == 0 then Right(None)
    else
      val denominator = ExactSum.zero()
      denominator.add(count.toDouble)
      val length = plan.channels * plan.channels * plan.binCount
      val outputReal = DVecBuilder.zeros(length)
      val outputImaginary = DVecBuilder.zeros(length)
      var index = 0
      var row = 0
      while row < plan.channels do
        var column = row
        while column < plan.channels do
          var k = 0
          while k < plan.binCount do
            val re = real(index).ratio(denominator).fold(e => throw new IllegalStateException(e.toString), identity)
            val im = imaginary(index).ratio(denominator).fold(e => throw new IllegalStateException(e.toString), identity)
            if !re.isFinite || !im.isFinite then return Left(BoundedCsdPlan.failure("nonfinite mean cross spectrum"))
            val upper = (row * plan.channels + column) * plan.binCount + k
            val lower = (column * plan.channels + row) * plan.binCount + k
            outputReal(upper) = re; outputReal(lower) = re
            outputImaginary(upper) = im; outputImaginary(lower) = if row == column then 0.0 else -im
            index += 1; k += 1
          column += 1
        row += 1
      FrequencyAxis.realFft(plan.stft.nfft, plan.stft.sampleRate).map(axis =>
        Some(new CsdMatrix(plan.channels, plan.binCount, outputReal.result(), outputImaginary.result(), count, axis)))
