package signal4s.multirate

import gale.linalg.Vec
import signal4s.*

class ResamplingBoundarySuite extends munit.FunSuite:
  test("output lengths use wide products and reject unrepresentable or invalid shapes"):
    assertEquals(Upfirdn.outputLength(1, 2, 1, Int.MaxValue), 1)
    assertEquals(Upfirdn.outputLength(3, 3, Int.MaxValue - 1, Int.MaxValue), 3)
    intercept[IllegalArgumentException](Upfirdn.outputLength(1, 2, Int.MaxValue, 1))
    intercept[IllegalArgumentException](Upfirdn.outputLength(1, 0, 1, 0))
    intercept[IllegalArgumentException](Upfirdn.outputLength(1, -1, 1, 1))
    assert(Upfirdn(Vec(1.0), Vec(1.0, 2.0), Int.MaxValue, 1).isLeft)

  test("tiny input and extreme factors preserve kept phases without giant phase banks"):
    assertEquals(Upfirdn(Vec(1.0), Vec(1.0, 2.0), 1, Int.MaxValue).orThrow.toSeq, Seq(1.0))
    assertEquals(Upfirdn(Vec(1.0, 2.0, 1.0), Vec(1.0, 2.0, 3.0), Int.MaxValue - 1, Int.MaxValue).orThrow.toSeq,
      Seq(1.0, 4.0, 3.0))

  test("delay metadata is in input-sample units on the reduced grid"):
    val prototype = Vec(0.25, 0.5, 0.25)
    assertEqualsDouble(PolyphaseResampler(prototype, 2, 3).orThrow.delay, 0.5, 0.0)
    assertEqualsDouble(PolyphaseResampler(prototype, 4, 6).orThrow.delay, 0.5, 0.0)
    assertEqualsDouble(PolyphaseResampler(prototype, 1, 2).orThrow.delay, 1.0, 0.0)

  test("wide counters retain ordinary partition, empty and flush lifecycle behavior"):
    val runner = PolyphaseResampler(Vec(1.0), 2, 3).orThrow
    assertEquals(runner.consume(Vec(1.0)).orThrow.toSeq, Seq(1.0))
    assertEquals(runner.samplesConsumed, 1L)
    assertEquals(runner.phase, 1)
    assertEquals(runner.consume(Vec.zeros(0)).orThrow.length, 0)
    assertEquals(runner.consume(Vec(2.0, 3.0)).orThrow.toSeq, Seq(0.0))
    assertEquals(runner.samplesConsumed, 3L)
    assertEquals(runner.flush().orThrow.length, 0)
    assert(runner.consume(Vec(4.0)).isLeft)
