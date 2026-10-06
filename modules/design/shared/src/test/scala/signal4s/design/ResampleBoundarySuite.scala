package signal4s.design

import gale.linalg.Vec
import signal4s.*

class ResampleBoundarySuite extends munit.FunSuite:
  test("valid empty designed resampling terminates across reduced ratios"):
    for (up, down) <- List((2, 3), (3, 2), (1, 5), (4, 6), (1, 1)) do
      assertEquals(ResamplePoly(Vec.zeros(0), up, down).orThrow.length, 0)
    assert(ResamplePoly(Vec.zeros(0), 0, 1).isLeft)

  test("unrepresentable resampling plans fail before prototype allocation"):
    assert(ResamplePoly(Vec(1.0, 2.0), Int.MaxValue, 1).isLeft)
    assert(ResamplePoly(Vec(1.0, 2.0), 1, Int.MaxValue).isLeft)

  test("one-sample and final short observations retain centered output length"):
    for (up, down) <- List((2, 3), (3, 2), (4, 6)) do
      val result = ResamplePoly(Vec(1.0), up, down).orThrow
      assertEquals(result.length, (up + down - 1) / down)
      assert(result.toSeq.forall(_.isFinite))
