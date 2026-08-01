package signal4s.laws

import gale.linalg.Vec
import signal4s.*
import signal4s.fft.*

class FftLawSuite extends munit.FunSuite:

  private def realSignal(n: Int): gale.linalg.DVec =
    Vec.tabulate(n)(i => math.sin(0.3 * i) + 0.1 * math.cos(0.7 * i))

  private def complexSignal(n: Int): ComplexVector =
    ComplexVector
      .tabulate(n)(i => Complex(math.sin(0.2 * i), math.cos(0.4 * i)))
      .orThrow

  test("complex round-trip power-of-two all norms"):
    List(
      FftNormalization.Backward,
      FftNormalization.Forward,
      FftNormalization.Orthonormal
    ).foreach { norm =>
      val plan = FftPlan(64, norm).orThrow
      FftLaws.complexRoundTrip(plan, complexSignal(64))
      FftLaws.parseval(plan, complexSignal(64))
    }

  test("complex round-trip mixed-radix 12 and 15"):
    val plan12 = FftPlan(12, FftNormalization.Backward).orThrow
    val plan15 = FftPlan(15, FftNormalization.Backward).orThrow
    FftLaws.complexRoundTrip(plan12, complexSignal(12))
    FftLaws.complexRoundTrip(plan15, complexSignal(15))

  test("complex round-trip prime Bluestein 17"):
    val plan = FftPlan(17, FftNormalization.Backward).orThrow
    FftLaws.complexRoundTrip(plan, complexSignal(17))
    FftLaws.parseval(plan, complexSignal(17))

  test("real round-trip and conjugate symmetry"):
    val plan = RealFftPlan(32, FftNormalization.Backward).orThrow
    val x = realSignal(32)
    FftLaws.realRoundTrip(plan, x)
    FftLaws.conjugateSymmetry(plan, x)

  test("real round-trip Bluestein length 17"):
    val plan = RealFftPlan(17, FftNormalization.Orthonormal).orThrow
    FftLaws.realRoundTrip(plan, realSignal(17))

  test("workspace reuse does not mutate prior ComplexVector results"):
    val plan = FftPlan(16, FftNormalization.Backward).orThrow
    val x = complexSignal(16)
    val ws = plan.newWorkspace()
    plan.forwardInto(x, ws).orThrow
    val y1 = plan.resultFrom(ws).orThrow
    plan.forwardInto(complexSignal(16).scale(2.0), ws).orThrow
    FftLaws.assertClose(y1, plan.forward(x).orThrow)
