package signal4s.fft

import gale.linalg.DVec
import signal4s.SampleRate
import munit.FunSuite

class MultitaperSuite extends FunSuite:
  private val budget=BigInt(268435456)
  private def ok[A](value: Either[?,A]): A=value.fold(e=>fail(e.toString),identity)
  private def tapers(c: MultitaperReferenceCase): DpssTapers=ok(DpssTapers.withinBudget(c.length,c.nw,c.candidates,c.floor,budget))
  private def plan(c: MultitaperReferenceCase): BoundedMultitaperPlan=
    ok(BoundedMultitaperPlan.withinBudget(2,tapers(c),ok(SampleRate.hertz(c.rate)),c.nfft,
      if c.detrend=="Constant" then Detrend.Mean else Detrend.None,budget))
  private def frame(c: MultitaperReferenceCase,index: Int): Vector[DVec]=
    Vector.tabulate(2)(row=>DVec.tabulate(c.length)(t=>c.inputs((index*2+row)*c.length+t)))

  for c<-MultitaperReference.cases do
    test(s"DPSS ${c.length}/${c.nw} matches pinned SciPy concentration, subspace and normalization"):
      val selected=tapers(c)
      assertEquals(selected.selectedIndices,c.selected)
      selected.candidateConcentrations.zip(c.concentrations).foreach((a,b)=>assertEqualsDouble(a,b,1e-10))
      for i<-selected.taps.indices do
        val expected=c.tapers.slice(c.selected(i)*c.length,(c.selected(i)+1)*c.length)
        val dot=expected.indices.map(k=>expected(k)*selected.taps(i)(k)).sum
        val sign=if dot<0 then -1.0 else 1.0
        expected.indices.foreach(k=>assertEqualsDouble(sign*selected.taps(i)(k),expected(k),1e-9))
        for j<-0 to i do
          val inner=(0 until c.length).map(k=>selected.taps(i)(k)*selected.taps(j)(k)).sum
          assertEqualsDouble(inner,if i==j then 1.0 else 0.0,1e-9)

    test(s"multichannel ${c.length}/${c.nfft} density matches independent direct FFT and Hermitian laws"):
      val p=plan(c);val state=ok(p.newAccumulator(budget))
      (0 until 3).foreach(i=>ok(state.accumulateFrame(frame(c,i))))
      val value=ok(state.finishOwned()).get
      assertEquals(value.acceptedWindows,3L)
      for i<-c.real.indices do
        assertEqualsDouble(value.real(i),c.real(i),1e-12+1e-9*math.abs(c.real(i)))
        assertEqualsDouble(value.imaginary(i),c.imaginary(i),1e-12+1e-9*math.abs(c.imaginary(i)))
      for k<-0 until value.binCount do
        assertEquals(value(0,1,k).real,value(1,0,k).real)
        assertEquals(value(0,1,k).imaginary,-value(1,0,k).imaginary)
        assertEquals(value(0,0,k).imaginary,0.0)
        val cross=value(0,1,k)
        assert(cross.real*cross.real+cross.imaginary*cross.imaginary<=value(0,0,k).real*value(1,1,k).real+1e-12)

  test("same-plan merges and frame permutations are exact; snapshots remain owned"):
    val c=MultitaperReference.cases.head;val p=plan(c)
    val serial=ok(p.newAccumulator(budget));val left=ok(p.newAccumulator(budget));val right=ok(p.newAccumulator(budget))
    for i<-0 until 3 do ok(serial.accumulateFrame(frame(c,i)))
    ok(left.accumulateFrame(frame(c,2)));ok(right.accumulateFrame(frame(c,1)));ok(right.accumulateFrame(frame(c,0)))
    ok(left.merge(right))
    val value=ok(left.finishOwned()).get;val expected=ok(serial.finishOwned()).get
    assertEquals(value.real.toSeq,expected.real.toSeq);assertEquals(value.imaginary.toSeq,expected.imaginary.toSeq)
    val snapshot=value.real.toSeq
    ok(left.accumulateFrame(frame(c,0)))
    assertEquals(value.real.toSeq,snapshot)
    assert(left.merge(left).isLeft);assert(left.merge(ok(plan(c).newAccumulator(budget))).isLeft)

  test("Parseval integrates the mean tapered energy with density units"):
    val c=MultitaperReference.cases.head;val p=plan(c);val input=frame(c,0)
    val state=ok(p.newAccumulator(budget));ok(state.accumulateFrame(input));val value=ok(state.finishOwned()).get
    for row<-0 until 2 do
      val energy=p.tapers.taps.map: taper =>
        val numerator=(0 until c.length).map(t=>math.pow(input(row)(t)*taper(t),2)).sum
        numerator/(0 until c.length).map(t=>taper(t)*taper(t)).sum
      val expected=energy.sum/energy.size
      val integrated=(0 until value.binCount).map(k=>value(row,row,k).real).sum*c.rate/c.nfft
      assertEqualsDouble(integrated,expected,1e-10)

  test("a later taper failure cannot commit a partial frame"):
    val c=MultitaperReference.cases.head;val ts=tapers(c);val rate=ok(SampleRate.hertz(c.rate))
    val p=ok(BoundedMultitaperPlan.withinBudget(2,ts,rate,c.nfft,Detrend.None,budget));val state=ok(p.newAccumulator(budget))
    ok(state.accumulateFrame(frame(c,0)));val before=ok(state.finishOwned()).get
    val amplitude=math.sqrt(Double.MaxValue/4)*math.sqrt(c.rate)/math.abs(ts.taps.head(0))
    val impulse=DVec.tabulate(c.length)(i=>if i==0 then amplitude else 0.0)
    val rows=Vector(impulse,impulse)
    val first=ok(BoundedCsdPlan.withinBudget(2,ok(Window.fromTaps(ts.taps.head,WindowConvention.Symmetric)),rate,c.nfft,budget)).newAccumulator(budget)
    assert(ok(first).accumulateFrame(rows).isRight,"first taper must succeed so this exercises a later failure")
    assert(state.accumulateFrame(rows).isLeft)
    assertEquals(state.acceptedWindows,1L)
    val after=ok(state.finishOwned()).get
    assertEquals(after.real.toSeq,before.real.toSeq);assertEquals(after.imaginary.toSeq,before.imaginary.toSeq)
    ok(state.accumulateFrame(frame(c,1)))
    assertEquals(state.acceptedWindows,2L)

  test("invalid requests and exact budgets refuse before construction or commit"):
    for (n,nw,k,floor)<-Vector((1,0.2,1,0.0),(32,Double.NaN,3,0.0),(32,16.0,3,0.0),
        (32,2.0,0,0.0),(32,2.0,33,0.0),(32,2.0,3,-0.1),(513,3.0,5,0.0)) do
      assert(DpssTapers.withinBudget(n,nw,k,floor,budget).isLeft)
    val c=MultitaperReference.cases.head;val bytes=ok(DpssTapers.inspect(c.length,c.nw,c.candidates,c.floor))
    assert(DpssTapers.withinBudget(c.length,c.nw,c.candidates,c.floor,bytes-1).isLeft)
    assert(DpssTapers.withinBudget(c.length,c.nw,c.candidates,1.0,budget).isLeft)
    val p=plan(c)
    assert(p.newAccumulator(p.resources.runtimeBytes-1).isLeft)
    val empty=ok(p.newAccumulator(p.resources.runtimeBytes));assertEquals(ok(empty.finishOwned()),None)
    assert(empty.accumulateFrame(Vector.empty).isLeft)
    assert(empty.accumulateFrame(Vector(DVec.tabulate(c.length)(_=>Double.NaN),frame(c,0)(1))).isLeft)
    assertEquals(empty.acceptedWindows,0L)

  test("zero complete frames produce valid zero density, not absent observations"):
    val c=MultitaperReference.cases.head;val p=plan(c);val state=ok(p.newAccumulator(budget))
    ok(state.accumulateFrame(Vector.fill(2)(DVec.zeros(c.length))))
    val result=ok(state.finishOwned()).get
    assertEquals(result.acceptedWindows,1L)
    assert(result.real.toSeq.forall(_==0));assert(result.imaginary.toSeq.forall(_==0))
