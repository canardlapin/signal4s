package signal4s.design
import gale.linalg.Vec
class ResamplePlanSuite extends munit.FunSuite:
  test("bounded arbitrary global-phase windows agree with centered full batches"):
    for (up,down)<-List((2,3),(3,2),(4,6),(1,5),(1,1)) do
      val p=ResamplePolyPlan(up,down).toOption.get
      val x=Vec.tabulate(173)(i=>math.sin(i*0.71)+0.1*i)
      val full=ResamplePoly(x,up,down).toOption.get
      for (start,count)<-List((0,1),(1,3),(7,11),(full.length-3,3),(full.length,0)) do
        val support=p.inputWindow(x.length,start,count).toOption.get
        val chunk=x.slice(support.start.toInt,support.start.toInt+support.count).copy
        val y=p.processWindow(chunk,support.start,x.length,start,count).toOption.get
        assertEquals(y.toSeq,full.slice(start,start+count).toSeq)
  test("Long phase geometry does not convert clocks to Double or overflow"):
    val p=ResamplePolyPlan(2,3).toOption.get
    val source=9007199254741201L
    val q=6004799503160701L
    val s=p.inputWindow(source,q,2).toOption.get
    assert(s.start>9007199254740000L && s.count<100)
    val x=Vec.tabulate(s.count)(_=>1.0)
    val y=p.processWindow(x,s.start,source,q,2).toOption.get
    assert(y.toSeq.forall(_.isFinite))
    assert(ResamplePolyPlan(3,2).toOption.get.outputLength(Long.MaxValue).isLeft)
  test("support and output capacities refuse before allocation"):
    val p=ResamplePolyPlan(1,3).toOption.get
    assert(p.inputWindow(Long.MaxValue,0L,Int.MaxValue).isLeft)
    assert(p.inputWindow(10L,0L,-1).isLeft)
    assert(p.inputWindow(10L,Long.MaxValue,1).isLeft)
    assertEquals(p.coefficientBytes,BigInt(61)*8)
    assertEquals(p.processAdditionalBytes(7).toOption.get,BigInt(112))
    assert(p.processAdditionalBytes(-1).isLeft)
  test("nonfinite input and late overflow expose no reusable partial output"):
    assert(ResamplePoly(Vec(Double.NaN),2,3).isLeft)
    assert(ResamplePoly(Vec(Double.PositiveInfinity),1,1).isLeft)
    val p=ResamplePolyPlan(2,1).toOption.get
    val x=Vec.tabulate(100)(_=>Double.MaxValue)
    assert(p.processWindow(x,0,100,0,200).isLeft)
    assertEquals(x(0),Double.MaxValue)
  test("mismatched support refuses and owned results survive later calls"):
    val p=ResamplePolyPlan(2,3).toOption.get
    val x=Vec.tabulate(100)(_=>1.0)
    val s=p.inputWindow(100,5,3).toOption.get
    assert(p.processWindow(x,0,100,5,3).isLeft)
    val data=x.slice(s.start.toInt,s.start.toInt+s.count).copy
    val y=p.processWindow(data,s.start,100,5,3).toOption.get
    val old=y.toSeq
    p.processWindow(Vec.zeros(s.count),s.start,100,5,3).toOption.get
    assertEquals(y.toSeq,old)
  test("source exclusions propagate conservative centered support including zero taps"):
    val p=ResamplePolyPlan(2,3).toOption.get
    assertEquals(p.affectedOutputs(100,40,41).toOption.get,(17L,37L))
    assert(p.affectedOutputs(100,50,50).isLeft)
    assertEquals(p.inputWindow(0,0,0).toOption.get.count,0)
