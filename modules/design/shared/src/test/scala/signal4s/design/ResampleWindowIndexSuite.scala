package signal4s.design

import gale.linalg.{DVec,Vec}

class ResampleWindowIndexSuite extends munit.FunSuite:
  private def ok[E,A](v: Either[E,A]): A=v.fold(e=>fail(e.toString),identity)
  private def bits(v: Double): Long=java.lang.Double.doubleToRawLongBits(v)
  // Independent contributor scan with unbounded absolute arithmetic; no support quotient/tap recurrence.
  private def oracle(plan: ResamplePolyPlan,input: DVec,inputStart: Long,outputStart: Long,count: Int): Vector[Long]=
    Vector.tabulate(count): q =>
      val clock=(BigInt(outputStart)+q)*plan.ratio.down+plan.halfLength
      var sum=0.0;var n=0
      while n<input.length do
        val tap=clock-(BigInt(inputStart)+n)*plan.ratio.up
        if tap>=0 && tap<plan.coefficients.length && input(n)!=0.0 then sum+=input(n)*plan.coefficients(tap.toInt)
        n+=1
      bits(sum)
  private def check(plan: ResamplePolyPlan,length: Long,start: Long,count: Int): Unit=
    val window=ok(plan.inputWindow(length,start,count))
    val input=Vec.tabulate(window.count)(i=>if i%9==0 then -0.0 else (i%19-9)*0.125)
    val before=input.toSeq.toVector
    val result=ok(plan.processWindow(input,window.start,length,start,count))
    assertEquals(result.toSeq.map(bits).toVector,oracle(plan,input,window.start,start,count))
    assertEquals(input.toSeq.toVector,before)
  test("ordinary and boundary windows preserve exact ascending contributor arithmetic across ratios"):
    for (up,down)<-Vector((1,1),(2,3),(3,2),(7,11),(13,4),(5,17));length<-Vector(1L,31L,127L,256L) do
      val plan=ok(ResamplePolyPlan(up,down));val total=ok(plan.outputLength(length))
      for start<-Vector(0L,total/2,math.max(0L,total-17)).distinct do check(plan,length,start,math.min(17L,total-start).toInt)
      check(plan,length,total,0)
  test("large absolute clocks retain exact support and explicit overflow fallback"):
    for (up,down)<-Vector((1,1),(1,2),(2,3),(7,11),(3,5)) do
      val plan=ok(ResamplePolyPlan(up,down));val total=ok(plan.outputLength(Long.MaxValue))
      for delta<-Vector(1L,20L,100L) do
        val start=total-delta;check(plan,Long.MaxValue,start,math.min(delta,7L).toInt)
