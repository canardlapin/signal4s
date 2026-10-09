package signal4s.design

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*
import signal4s.fft.{WindowConvention, WindowSpec}
import signal4s.multirate.RateRatio

/** A bounded contiguous input range in the original segment's global ordinals. */
final case class ResampleInputWindow(start: Long, count: Int)

/** Immutable centered zero-padded rational resampling. The same designed prototype
  * and absolute output phase serve batch and arbitrary window evaluation.
  * Numeric capacity descriptions exclude caller input, objects/BigInt/allocator/GC.
  */
final class ResamplePolyPlan private (val ratio: RateRatio, val coefficients: DVec):
  val halfLength: Int = (coefficients.length-1)/2
  val coefficientBytes: BigInt = BigInt(coefficients.length)*8
  /** Conservative builder plus completed owned output payload; no source-length allocation. */
  def processAdditionalBytes(count: Int): Either[SignalError,BigInt] =
    if count<0 then fail("negative output capacity") else Right(BigInt(count)*16)
  val inputHalo: Int = ((halfLength.toLong+ratio.up-1)/ratio.up).toInt
  def maximumInputCount(outputCount: Int): Either[SignalError,Int] =
    if outputCount<0 then fail("negative output count")
    else if outputCount==0 then Right(0)
    else
      val n=((BigInt(outputCount-1)*ratio.down+2L*halfLength+ratio.up-1)/ratio.up)+1
      if n>Int.MaxValue then fail("input support exceeds bounded Int capacity") else Right(n.toInt)
  def outputLength(inputLength: Long): Either[SignalError,Long] =
    if inputLength<0 then fail("negative input length")
    else
      val n=(BigInt(inputLength)*ratio.up+ratio.down-1)/ratio.down
      if !n.isValidLong then fail("output length exceeds Long capacity") else Right(n.toLong)
  private def fail[A](message: String): Either[SignalError,A] =
    Left(SignalError.NumericalFailure("ResamplePolyPlan",message))
  private def floor(n: BigInt,d: Int): BigInt =
    val q=n/d
    if n.signum<0 && n%d!=0 then q-1 else q
  private def ceil(n: BigInt,d: Int): BigInt = -floor(-n,d)

  /** All observed contributors, including zero taps. Before/after segment is zero padding. */
  def inputWindow(inputLength: Long, outputStart: Long, outputCount: Int): Either[SignalError,ResampleInputWindow] =
    outputLength(inputLength).flatMap: total =>
      if outputStart<0 || outputCount<0 || BigInt(outputStart)+outputCount>total then fail("output window outside segment")
      else if outputCount==0 then Right(ResampleInputWindow(0L,0))
      else
        val lo=ceil(BigInt(outputStart)*ratio.down-halfLength,ratio.up).max(0).min(inputLength)
        val hi=(floor((BigInt(outputStart)+outputCount-1)*ratio.down+halfLength,ratio.up)+1).max(lo).min(inputLength)
        val count=hi-lo
        if count>Int.MaxValue then fail("input support exceeds bounded Int capacity")
        else Right(ResampleInputWindow(lo.toLong,count.toInt))

  /** Conservative affected output interval for a source exclusion [start,stop). */
  def affectedOutputs(inputLength: Long,start: Long,stop: Long): Either[SignalError,(Long,Long)] =
    outputLength(inputLength).flatMap: total =>
      if start<0 || stop<=start || stop>inputLength then fail("invalid source exclusion")
      else
        val a=ceil(BigInt(start)*ratio.up-halfLength,ratio.down).max(0).min(total)
        val b=(floor(BigInt(stop-1)*ratio.up+halfLength,ratio.down)+1).max(a).min(total)
        Right((a.toLong,b.toLong))

  /** Input must be exactly inputWindow's completed owned support, not a local phase reset.
    * Finite admission and late overflow return no partial reusable output. Ascending
    * input accumulation matches the batch upfirdn scatter order under ordinary rounding.
    */
  def processWindow(input: DVec,inputStart: Long,inputLength: Long,outputStart: Long,outputCount: Int): Either[SignalError,DVec] =
    inputWindow(inputLength,outputStart,outputCount).flatMap: required =>
      if inputStart!=required.start || input.length!=required.count then fail("completed input differs from required global support")
      else
        var i=0
        var finite=true
        while i<input.length && finite do {finite=input(i).isFinite;i+=1}
        if !finite then fail("nonfinite input")
        else
          val out=DVecBuilder.zeros(outputCount)
          var q=0
          var error=false
          while q<outputCount && !error do
            val clock=(BigInt(outputStart)+q)*ratio.down+halfLength
            val first=ceil(clock-(coefficients.length-1),ratio.up).max(required.start)
            val stop=(floor(clock,ratio.up)+1).min(BigInt(required.start)+required.count)
            var n=(first-required.start).toInt
            val end=(stop-required.start).toInt
            var tap=(clock-first*ratio.up).toInt
            var sum=0.0
            while n<end && !error do
              val value=input(n)
              // Match upfirdn's zero-input shortcut, including signed zero behavior.
              if value!=0.0 then sum+=value*coefficients(tap)
              error = !sum.isFinite
              n+=1
              tap-=ratio.up
            if !error then out(q)=sum
            q+=1
          if error then fail("nonfinite output/intermediate") else Right(out.result())

object ResamplePolyPlan:
  def apply(up: Int,down: Int,window: WindowSpec=WindowSpec.Kaiser(1,5.0,WindowConvention.Symmetric)): Either[SignalError,ResamplePolyPlan] =
    RateRatio(up,down).flatMap: ratio =>
      if ratio.isIdentity then Right(new ResamplePolyPlan(ratio,DVec.tabulate(1)(_=>1.0)))
      else
        val maxRate=math.max(ratio.up,ratio.down)
        val length=BigInt(20)*maxRate+1
        if length>Int.MaxValue then Left(SignalError.NumericalFailure("ResamplePolyPlan","prototype length exceeds Int capacity"))
        else
          FirDesign.lowPass(length.toInt,Frequency.unsafe(1.0/maxRate),SampleRate.unsafe(2.0),window).flatMap: designed =>
            val taps=DVec.tabulate(length.toInt)(i=>designed.fir.taps(i)*ratio.up)
            if !taps.toSeq.forall(_.isFinite) then Left(SignalError.NumericalFailure("ResamplePolyPlan","nonfinite prototype"))
            else Right(new ResamplePolyPlan(ratio,taps))
