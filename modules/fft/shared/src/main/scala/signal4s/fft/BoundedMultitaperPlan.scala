package signal4s.fft

import gale.linalg.DVec
import gale.numeric.ExactSum
import signal4s.{SampleRate,SignalError}

final case class MultitaperResources(runtimeBytes: BigInt,taperConstructionBytes: BigInt):
  val peakBytes: BigInt=runtimeBytes.max(taperConstructionBytes)

/** Fixed equal-weight taper means, then arithmetic common-window means.
  * One-sided density, including complex conj(first)*second CSD. Adaptive
  * weights, jackknife intervals and an inferential degrees of freedom are absent.
  */
final class BoundedMultitaperPlan private (val channels: Int,val tapers: DpssTapers,
    val nfft: Int,val sampleRate: SampleRate,val detrend: Detrend,val resources: MultitaperResources,
    private[fft] val frames: Vector[BoundedStftPlan],private[fft] val aggregation: BoundedCsdPlan):
  val binCount: Int=nfft/2+1
  val weighting: String="EqualWithinFrameThenEqualWindows"
  private[fft] val cells: Int=aggregation.cells
  def newAccumulator(numericBudgetBytes: BigInt): Either[SignalError,MultitaperAccumulator] =
    for
      _ <- if numericBudgetBytes>=resources.runtimeBytes then Right(()) else Left(BoundedCsdPlan.failure("multitaper runtime exceeds budget"))
      aggregate <- aggregation.newAccumulator(numericBudgetBytes)
      workspaces <- frames.foldLeft[Either[SignalError,Vector[StftWorkspace]]](Right(Vector.empty)): (state,plan) =>
        for previous<-state;workspace<-plan.newWorkspace(numericBudgetBytes) yield previous:+workspace
    yield new MultitaperAccumulator(this,aggregate,workspaces)

object BoundedMultitaperPlan:
  def inspect(channels: Int,tapers: DpssTapers,nfft: Int): Either[SignalError,MultitaperResources] =
    for
      base <- BoundedCsdPlan.inspect(channels,tapers.length,nfft)
      frame <- BoundedStftPlan.inspect(tapers.length,nfft)
      cells=BigInt(channels)*(channels+1)/2*(nfft/2+1)
      staging=cells*(tapers.count+1)*16+BigInt(channels)*(nfft/2+1)*32+
        ExactSum.resources.stateBytes*3+ExactSum.resources.ratioScratchBytes
    yield MultitaperResources(base.peakBytes+tapers.retainedBytes+staging+
      tapers.count*frame.peakBytes.max(frame.constructionAdditionalBytes+frame.windowBytes),tapers.constructionBytes)

  def withinBudget(channels: Int,tapers: DpssTapers,sampleRate: SampleRate,nfft: Int,
      detrend: Detrend,numericBudgetBytes: BigInt): Either[SignalError,BoundedMultitaperPlan] =
    for
      resource <- inspect(channels,tapers,nfft)
      _ <- if resource.peakBytes<=numericBudgetBytes then Right(()) else Left(BoundedCsdPlan.failure("multitaper construction/runtime exceeds budget"))
      _ <- if detrend==Detrend.None || detrend==Detrend.Mean then Right(()) else Left(BoundedCsdPlan.failure("multitaper v1 supports None or Mean (constant) detrending"))
      windows <- tapers.taps.foldLeft[Either[SignalError,Vector[Window]]](Right(Vector.empty)): (state,taps) =>
        for previous<-state;window<-Window.fromTaps(taps,WindowConvention.Symmetric) yield previous:+window
      frames <- windows.foldLeft[Either[SignalError,Vector[BoundedStftPlan]]](Right(Vector.empty)): (state,window) =>
        for previous<-state;plan<-BoundedStftPlan.withinBudget(window,sampleRate,nfft,numericBudgetBytes,detrend,StftScaling.Density) yield previous:+plan
      aggregate <- BoundedCsdPlan.withinBudget(channels,windows.head,sampleRate,nfft,numericBudgetBytes,detrend)
    yield new BoundedMultitaperPlan(channels,tapers,nfft,sampleRate,detrend,resource,frames,aggregate)

/** Single-owner transaction: every taper must complete before the frame commits. */
final class MultitaperAccumulator private[fft] (val plan: BoundedMultitaperPlan,
    private val aggregate: CsdAccumulator,workspaces: Vector[StftWorkspace]) extends CsdFrameAccumulator:
  private val real=Array.fill(plan.tapers.count)(new Array[Double](plan.cells))
  private val imaginary=Array.fill(plan.tapers.count)(new Array[Double](plan.cells))
  private val meanReal=new Array[Double](plan.cells)
  private val meanImaginary=new Array[Double](plan.cells)
  def acceptedWindows: Long=aggregate.acceptedWindows

  def accumulateFrame(rows: IndexedSeq[DVec]): Either[SignalError,Unit] =
    if rows.size!=plan.channels || rows.exists(_.length!=plan.tapers.length) then
      return Left(BoundedCsdPlan.failure("multitaper requires a complete common frame in all rows"))
    if acceptedWindows==BoundedCsdPlan.MaxWindows then return Left(BoundedCsdPlan.failure("exact window count capacity exceeded"))
    var taper=0
    while taper<plan.tapers.count do
      val spectra=new Array[ComplexVector](plan.channels)
      var row=0
      while row<plan.channels do
        plan.frames(taper).analyzeInto(rows(row),workspaces(taper)) match
          case Left(error)=>return Left(error)
          case Right(frame)=>spectra(row)=frame.bins
        row+=1
      CrossPeriodogram.write(spectra,plan.nfft,real(taper),imaginary(taper)) match
        case Left(error)=>return Left(error)
        case Right(_)=>()
      taper+=1
    val denominator=ExactSum.zero();denominator.add(plan.tapers.count.toDouble)
    var cell=0
    while cell<plan.cells do
      val re=ExactSum.zero();val im=ExactSum.zero()
      taper=0
      while taper<plan.tapers.count do
        re.add(real(taper)(cell)).fold(e=>throw new IllegalStateException(e.message),identity)
        im.add(imaginary(taper)(cell)).fold(e=>throw new IllegalStateException(e.message),identity)
        taper+=1
      re.ratio(denominator) match
        case Left(e)=>return Left(BoundedCsdPlan.failure(e.message))
        case Right(value)=>meanReal(cell)=value
      im.ratio(denominator) match
        case Left(e)=>return Left(BoundedCsdPlan.failure(e.message))
        case Right(value)=>meanImaginary(cell)=value
      cell+=1
    aggregate.accumulateSpectrum(meanReal,meanImaginary)

  def merge(other: MultitaperAccumulator): Either[SignalError,Unit] =
    if this eq other then Left(BoundedCsdPlan.failure("self merge duplicates observations"))
    else if plan ne other.plan then Left(BoundedCsdPlan.failure("different multitaper plan identity"))
    else aggregate.merge(other.aggregate)

  def finishOwned(): Either[SignalError,Option[CsdMatrix]]=aggregate.finishOwned()
