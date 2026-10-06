package signal4s.multirate

import gale.linalg.Vec
import java.lang.management.ManagementFactory
import signal4s.*

class ResamplingClockAllocationSuite extends munit.FunSuite:
  private def seed(runner: PolyphaseResampler, name: String, value: Long): Unit =
    val field = runner.getClass.getDeclaredField(name)
    field.setAccessible(true)
    field.setLong(runner, value)

  test("long history crosses the Int boundary without corrupting decimation phase"):
    val runner = PolyphaseResampler(Vec(1.0), 1, 3).orThrow
    seed(runner, "inputCount", Int.MaxValue.toLong)
    seed(runner, "upIndex", Int.MaxValue.toLong)
    val output = runner.consume(Vec(1.0, 2.0, 3.0)).orThrow
    assertEquals(output.toSeq, Seq(3.0))
    assertEquals(runner.samplesConsumed, Int.MaxValue.toLong + 3)

  test("clock exhaustion refuses consume and flush before changing state"):
    val runner = PolyphaseResampler(Vec(1.0), 1, 3).orThrow
    seed(runner, "inputCount", Long.MaxValue)
    seed(runner, "upIndex", Long.MaxValue)
    assert(runner.consume(Vec(9.0)).isLeft)
    assertEquals(runner.samplesConsumed, Long.MaxValue)
    assert(!runner.isFlushed)
    assertEquals(runner.consume(Vec.zeros(0)).orThrow.length, 0)
    val draining = PolyphaseResampler(Vec(1.0, 1.0), 1, 3).orThrow
    seed(draining, "inputCount", Long.MaxValue)
    seed(draining, "upIndex", Long.MaxValue)
    assert(draining.flush().isLeft)
    assert(!draining.isFlushed)

  test("tiny decimation allocation stays bounded independently of an unused phase denominator"):
    val available = ManagementFactory.getThreadMXBean match
      case value: com.sun.management.ThreadMXBean if value.isThreadAllocatedMemorySupported => Some(value)
      case _ => None
    assume(available.nonEmpty, "thread allocation measurement unavailable")
    val bean = available.get
    val enabled = bean.isThreadAllocatedMemoryEnabled
    if !enabled then bean.setThreadAllocatedMemoryEnabled(true)
    try
      val h = Vec(0.25, 0.5, 0.25)
      val x = Vec(1.0, 2.0)
      var checksum = 0.0
      for _ <- 0 until 256 do checksum += Upfirdn(h, x, 1, 4096).orThrow(0)
      val thread = Thread.currentThread().threadId()
      val before = bean.getThreadAllocatedBytes(thread)
      for _ <- 0 until 128 do checksum += Upfirdn(h, x, 1, 4096).orThrow(0)
      val allocated = bean.getThreadAllocatedBytes(thread) - before
      assertEqualsDouble(checksum, 96.0, 0.0)
      assert(allocated >= 0 && allocated < 128L * 8192,
        s"$allocated bytes for128 tiny calls; allocating phase banks proportional to down is not bounded by the input")
    finally
      if !enabled then bean.setThreadAllocatedMemoryEnabled(false)
