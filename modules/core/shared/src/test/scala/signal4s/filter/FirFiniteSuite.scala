package signal4s.filter

import gale.linalg.{Vec,MutableDVec}
import signal4s.*

class FirFiniteSuite extends munit.FunSuite:
  test("FIR refuses nonfinite coefficients"):
    assert(Fir.causal(Vec(1,Double.NaN)).isLeft)

  test("late overflowing output refuses without changing state or destination"):
    val r=Fir.causal(Vec(2,1)).orThrow.newRunner()
    r.process(Vec(1)).orThrow
    val before=r.snapshot.values.toSeq
    val out=MutableDVec.zeros(2)
    out(0)=123;out(1)=123
    assert(r.processInto(Vec(1,Double.MaxValue),out).isLeft)
    assertEquals(r.snapshot.values.toSeq,before)
    assertEquals(out.toVec.toSeq,Seq(123.0,123.0))

  test("nonfinite input and state restore preserve a completed prefix"):
    val f=Fir.causal(Vec(0.25,0.5,0.25)).orThrow
    val r=f.newRunner();r.process(Vec(1,2)).orThrow
    val before=r.snapshot.values.toSeq
    for bad<-Vector(Double.NaN,Double.PositiveInfinity,Double.NegativeInfinity) do
      assert(r.process(Vec(1,bad)).isLeft)
      assert(r.restore(FilterState.from(Vec(1,bad))).isLeft)
      assert(f.newRunner(FilterState.from(Vec(bad,1))).isLeft)
      assertEquals(r.snapshot.values.toSeq,before)

  test("capacities match delay snapshots and staged/owned array payloads"):
    val f=Fir.causal(Vec(0.25,0.5,0.25)).orThrow
    assertEquals(f.resources.coefficientBytes,BigInt(32))
    assertEquals(f.resources.stateBytes,BigInt(f.newRunner().snapshot.length)*8)
    assertEquals(f.resources.snapshotBytes,BigInt(16))
    assertEquals(f.resources.processIntoScratchBytes(5),Right(BigInt(56)))
    assertEquals(f.resources.ownedProcessAdditionalBytes(5),Right(BigInt(136)))
    assert(f.resources.processIntoScratchBytes(-1).isLeft)

  test("partitioned ordinary finite outputs and owned snapshots retain causal arithmetic"):
    val f=Fir.causal(Vec(0.25,0.5,0.25)).orThrow
    val x=Vec.tabulate(17)(i=>math.sin(i*0.7)+0.1*i)
    val expected=f.process(x).orThrow
    val r=f.newRunner();val a=r.process(x.slice(0,3)).orThrow;val b=r.process(x.slice(3,11)).orThrow;val c=r.process(x.slice(11,17)).orThrow
    assertEquals(a.toSeq++b.toSeq++c.toSeq,expected.toSeq)
    r.reset();assertEquals(a.toSeq,expected.slice(0,3).toSeq)

  test("constructor refuses empty kernels and finite refusal retains input"):
    assert(Fir.causal(Vec.zeros(0)).isLeft)
    val x=Vec(1,Double.MaxValue);val before=x.toSeq
    assert(Fir.causal(Vec(2,1)).orThrow.process(x).isLeft)
    assertEquals(x.toSeq,before)
