package signal4s.design
import gale.linalg.Vec
class ResampleFiniteReproSuite extends munit.FunSuite:
  test("nonfinite resampling input is a typed refusal"):
    assert(ResamplePoly(Vec(Double.NaN),2,3).isLeft)
