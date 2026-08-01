package signal4s.fft

import gale.linalg.Vec
import signal4s.*

class FftSuite extends munit.FunSuite:

  test("delta spectrum is flat under Backward"):
    val plan = FftPlan(8, FftNormalization.Backward).orThrow
    val x = ComplexVector.tabulate(8)(i => if i == 0 then Complex.One else Complex.Zero).orThrow
    val X = plan.forward(x).orThrow
    var k = 0
    while k < 8 do
      assertEqualsDouble(X.real(k), 1.0, 1e-12)
      assertEqualsDouble(X.imaginary(k), 0.0, 1e-12)
      k += 1

  test("RealFftPlan bin count and axis"):
    val fs = SampleRate.hertz(1000.0).orThrow
    val plan = RealFftPlan(16, FftNormalization.Backward, fs).orThrow
    assertEquals(plan.binCount, 9)
    val x = Vec.tabulate(16)(_.toDouble)
    val spec = plan.forward(x).orThrow
    assertEquals(spec.bins.length, 9)
    assertEqualsDouble(spec.frequencies.frequencyAt(1).hertz, 1000.0 / 16.0, 1e-12)

  test("rejects non-positive length"):
    assert(FftPlan(0, FftNormalization.Backward).isLeft)
    assert(RealFftPlan(-1, FftNormalization.Forward).isLeft)
