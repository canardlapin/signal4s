package signal4s.multirate

import gale.linalg.Vec
import signal4s.*

class ResamplerCheckpointSuite extends munit.FunSuite:
  private val taps = Vec(0.25, 0.5, 0.25)

  test("owned checkpoints reproduce continuation and retain ring position and global phase"):
    for (up, down) <- List((2, 3), (3, 2), (5, 4)) do
      val runner = PolyphaseResampler(taps, up, down).orThrow
      val _ = runner.consume(Vec(1.0, 2.0, 3.0, 4.0, 5.0)).orThrow
      val state = runner.snapshot
      val registers = state.registers.toSeq
      val expected = runner.consume(Vec(6.0, 7.0, 8.0)).orThrow
      val expectedTail = runner.flush().orThrow
      val restored = PolyphaseResampler(taps.copy, up, down).orThrow
      restored.restore(state).orThrow
      assertEquals(restored.samplesConsumed, 5L)
      assertEquals(restored.consume(Vec(6.0, 7.0, 8.0)).orThrow.toSeq, expected.toSeq)
      assertEquals(restored.flush().orThrow.toSeq, expectedTail.toSeq)
      assertEquals(state.registers.toSeq, registers)

  test("checkpoint mismatch refuses restore before mutating an existing stream"):
    val source = PolyphaseResampler(taps, 2, 3).orThrow
    val _ = source.consume(Vec(1.0, 2.0)).orThrow
    val state = source.snapshot
    for runner <- List(PolyphaseResampler(taps, 3, 2).orThrow,
        PolyphaseResampler(Vec(0.5, 0.5, 0.0), 2, 3).orThrow) do
      val _ = runner.consume(Vec(9.0)).orThrow
      val before = runner.snapshot
      assert(runner.restore(state).isLeft)
      assertEquals(runner.samplesConsumed, before.samplesConsumed)
      assertEquals(runner.snapshot.registers.toSeq, before.registers.toSeq)
      assertEquals(runner.isFlushed, before.isFlushed)

  test("flushed snapshots remain closed, while reset explicitly starts a new stream"):
    val runner = PolyphaseResampler(taps, 2, 3).orThrow
    val first = runner.consume(Vec(1.0, 2.0)).orThrow
    val retained = first.toSeq
    val _ = runner.flush().orThrow
    val closed = runner.snapshot
    runner.reset()
    assertEquals(runner.samplesConsumed, 0L)
    assertEquals(runner.consume(Vec(1.0, 2.0)).orThrow.toSeq, retained)
    runner.restore(closed).orThrow
    assert(runner.consume(Vec(1.0)).isLeft)
    assertEquals(runner.flush().orThrow.length, 0)
