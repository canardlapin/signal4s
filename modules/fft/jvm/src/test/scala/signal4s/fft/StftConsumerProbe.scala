package signal4s.fft
import gale.linalg.Vec
import signal4s.*
object StftConsumerProbe:
 def main(args:Array[String]):Unit=
  val fs=SampleRate.hertz(16).orThrow
  val hann=Window.fromSpec(WindowSpec.Hann(4,WindowConvention.Periodic)).orThrow
  val rect=Window.fromSpec(WindowSpec.Rectangular(4,WindowConvention.Periodic)).orThrow
  val sparse=StftPlan(hann,4,fs,4)
  println(s"hann_hop4_constructor=$sparse")
  sparse.foreach: plan =>
   val input=Vec.tabulate(8)(i=>if i==2 then 1.0 else 0.0)
   val actual=plan.synthesize(plan.analyze(input).orThrow,8).orThrow
   println(s"hann_hop4_lost_impulse_input=${input(2)} output=${actual(2)}")
  println(s"rect_hop5_dual=${CanonicalDual.fromAnalysis(rect,5)}")
  val count=FramePlan(4,2,FrameAlignment.Centered).orThrow.frameCount(Int.MaxValue-2)
  println(s"centered_large_length_frame_count=$count")
  val fromStart=StftPlan(rect,2,fs,4,FrameAlignment.FromStart).orThrow
  println(s"from_start_first_time=${fromStart.analyze(Vec.tabulate(8)(_.toDouble)).orThrow.times.first.value}")
