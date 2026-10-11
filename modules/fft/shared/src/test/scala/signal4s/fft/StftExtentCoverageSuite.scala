package signal4s.fft
import gale.linalg.Vec
import signal4s.*

class StftExtentCoverageSuite extends munit.FunSuite:
 private val fs=SampleRate.hertz(16).orThrow
 private def rect(n:Int)=Window.fromSpec(WindowSpec.Rectangular(n,WindowConvention.Periodic)).orThrow

 test("zero overlap phases cannot acquire a canonical reconstruction dual"):
  val hann=Window.fromSpec(WindowSpec.Hann(4,WindowConvention.Periodic)).orThrow
  assert(CanonicalDual.fromAnalysis(hann,4).isLeft)
  assert(StftPlan(hann,4,fs,4).isLeft)
  assert(CanonicalDual.fromAnalysis(rect(4),5).isLeft)
  val zero=Window.fromTaps(Vec(0.0,0.0,0.0,0.0),WindowConvention.Periodic).orThrow
  assert(CanonicalDual.fromAnalysis(zero,2).isLeft)

 test("canonical dual is invariant to representable window scale"):
  for scale<-Vector(1e-280,1e-100,1e-8,1.0,1e100) do
   val window=Window.fromTaps(Vec.tabulate(4)(_=>scale),WindowConvention.Periodic).orThrow
   val dual=CanonicalDual.fromAnalysis(window,2).orThrow
   for i<-0 until 4 do assertEqualsDouble(dual.taps(i)*scale,0.5,1e-14)

 test("weak covered residues remain observable beside strong residues"):
  val window=Window.fromTaps(Vec(1e-140,1.0,1e-140,1.0),WindowConvention.Periodic).orThrow
  val dual=CanonicalDual.fromAnalysis(window,2).orThrow
  for residue<-0 until 2 do
   assertEqualsDouble(window.taps(residue)*dual.taps(residue)+window.taps(residue+2)*dual.taps(residue+2),1.0,1e-14)

 test("unrepresentable dual coefficients refuse rather than returning infinity"):
  val window=Window.fromTaps(Vec.tabulate(4)(_=>Double.MinPositiveValue),WindowConvention.Periodic).orThrow
  assert(CanonicalDual.fromAnalysis(window,2).isLeft)

 test("supported overlapping Hann reconstructs every impulse position"):
  val window=Window.fromSpec(WindowSpec.Hann(8,WindowConvention.Periodic)).orThrow
  val plan=StftPlan(window,4,fs,8).orThrow
  for position<-0 until 13 do
   val input=Vec.tabulate(13)(i=>if i==position then 1.0 else 0.00)
   val actual=plan.synthesize(plan.analyze(input).orThrow,13).orThrow
   for i<-0 until 13 do assertEqualsDouble(actual(i),input(i),1e-12)

 test("representable near-Int-limit counts remain positive without allocating input"):
  val frames=FramePlan(4,2,FrameAlignment.Centered).orThrow
  assertEquals(frames.frameCount(Int.MaxValue-2),1073741824)
  assertEquals(frames.frameCountChecked(Int.MaxValue.toLong-2).orThrow,1073741824L)

 test("Long and legacy Int extents refuse unrepresentable results explicitly"):
  val frames=FramePlan(4,1,FrameAlignment.Centered).orThrow
  assert(frames.frameCountChecked(Long.MaxValue).isLeft)
  assert(frames.frameCountChecked(-1).isLeft)
  intercept[ArithmeticException](frames.frameCount(Int.MaxValue))
  val spaced=FramePlan(4,3,FrameAlignment.FromStart).orThrow
  assert(spaced.frameStartChecked(Long.MaxValue).isLeft)
  intercept[ArithmeticException](spaced.frameStart(Int.MaxValue))

 test("time columns label even and odd discrete centers for both alignments"):
  for length<-Vector(3,4);alignment<-Vector(FrameAlignment.Centered,FrameAlignment.FromStart) do
   val plan=StftPlan(rect(length),1,fs,length,alignment).orThrow
   val tf=plan.analyze(Vec.tabulate(8)(_.toDouble)).orThrow
   val first=plan.frames.firstFrameStart+length/2
   for column<-0 until tf.frameCount do
    assertEqualsDouble(tf.times.timeAt(column).value,(first+column).toDouble/16,1e-15)
