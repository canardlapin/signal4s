package signal4s.fft

import gale.linalg.{DVec, Vec}
import signal4s.*

class BoundedMorletSuite extends munit.FunSuite:
  private val fs = SampleRate.hertz(128).orThrow
  private val budget = BigInt(64) * 1024 * 1024
  private def plan = BoundedMorletPlan.withinBudget(fs, Vector(8,16,32), Vector(3,5,7), budget).orThrow

  test("public MNE kernels and full-support impulse/tone coefficients preserve phase and amplitude"):
    for fixture <- MorletGoldens.cases do
      val p = BoundedMorletPlan.withinBudget(fs, Vector(fixture.frequency), Vector(fixture.cycles), budget, fixture.zeroMean).orThrow
      val kernel = p.kernels.head
      assertEquals(kernel.taps.length, fixture.kernelReal.size)
      for i <- 0 until kernel.taps.length do
        assert(math.abs(kernel.taps.real(i) - fixture.kernelReal(i)) <= 1e-14 + 1e-12 * math.abs(fixture.kernelReal(i)), fixture.name)
        assert(math.abs(kernel.taps.imaginary(i) - fixture.kernelImaginary(i)) <= 1e-14 + 1e-12 * math.abs(fixture.kernelImaginary(i)), fixture.name)
      val input = DVec.fromSeq(fixture.input)
      val result = p.analyzeBlock(input, 0, input.length, budget).orThrow
      val maximum = fixture.real.zip(fixture.imaginary).map((a,b) => math.hypot(a,b)).max
      for time <- 0 until input.length do
        val valid = time >= kernel.halfWidth && time + kernel.halfWidth < input.length
        assertEquals(result.isValid(0,time), valid)
        if valid then
          assert(math.abs(result.real(time) - fixture.real(time)) <= 1e-30 + 1e-10 * maximum, (fixture.name,time))
          assert(math.abs(result.imaginary(time) - fixture.imaginary(time)) <= 1e-30 + 1e-10 * maximum, (fixture.name,time))
        else
          assertEquals(result.real(time), 0.0)
          assertEquals(result.imaginary(time), 0.0)

  test("haloed chunking gives identical owned coefficients and frequency-specific edge masks"):
    val p = plan
    val input = Vec.tabulate(256)(i => math.sin(i * 0.1) * 1e-6)
    val whole = p.analyzeBlock(input, 0, input.length, budget).orThrow
    for start <- 0 until 256 by 17 do
      val count = math.min(17,256-start)
      val low = math.max(0,start-p.maximumHalo)
      val high = math.min(256,start+count+p.maximumHalo)
      val part = DVec.tabulate(high-low)(i => input(low+i))
      val result = p.analyzeBlock(part,start-low,count,budget).orThrow
      for f <- p.kernels.indices; i <- 0 until count do
        assertEquals(result.isValid(f,i),whole.isValid(f,start+i))
        assertEquals(result.real(f*count+i),whole.real(f*256+start+i))
        assertEquals(result.imaginary(f*count+i),whole.imaginary(f*256+start+i))
    val before = whole.real.toSeq
    p.analyzeBlock(Vec.fill(256)(0.0),0,256,budget).orThrow
    assertEquals(whole.real.toSeq,before)

  test("nonfinite input and unrepresentable convolution refuse without publishing partial output"):
    val p = plan
    assert(p.analyzeBlock(Vec.tabulate(256)(i => if i==255 then Double.NaN else 1.0),0,256,budget).isLeft)
    val dc = BoundedMorletPlan.withinBudget(fs,Vector(8),Vector(1),budget,zeroMean=false).orThrow
    assert(dc.analyzeBlock(Vec.fill(256)(Double.MaxValue),0,256,budget).isLeft)
    assert(p.analyzeBlock(Vec.fill(256)(1.0),0,256,p.peakBytes(256)-1).isLeft)
    assert(p.analyzeBlock(Vec.fill(256)(1.0),250,7,budget).isLeft)
    assert(p.analyzeBlock(Vec.fill(256)(1.0),0,-1,budget).isLeft)

  test("frequency cycles support and construction budgets are explicit"):
    for f <- Vector(0.0,-1.0,64.0,Double.NaN,Double.PositiveInfinity) do
      assert(BoundedMorletPlan.withinBudget(fs,Vector(f),Vector(3),budget).isLeft)
    for c <- Vector(0.0,0.5,65.0,Double.NaN) do
      assert(BoundedMorletPlan.withinBudget(fs,Vector(8),Vector(c),budget).isLeft)
    assert(BoundedMorletPlan.withinBudget(fs,Vector(16,8),Vector(3),budget).isLeft)
    assert(BoundedMorletPlan.withinBudget(fs,Vector(8,8),Vector(3),budget).isLeft)
    assert(BoundedMorletPlan.withinBudget(fs,Vector(1e-9),Vector(3),budget).isLeft)
    assert(BoundedMorletPlan.withinBudget(fs,Vector(8),Vector(3),0).isLeft)
    assert(BoundedMorletPlan.withinBudget(fs,Vector(8,16),Vector(3),budget).isRight)
