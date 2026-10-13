package signal4s.filter

import gale.linalg.{MutableDVec, Vec}
import java.lang.management.ManagementFactory
import signal4s.*

class FirAllocationSuite extends munit.FunSuite:
  test("warmed caller-output FIR allocates no per-sample payload"):
    val bean = ManagementFactory.getThreadMXBean match
      case value: com.sun.management.ThreadMXBean if value.isThreadAllocatedMemorySupported => value
      case _ => fail("HotSpot allocation counter required for this JVM regression")
    val enabled = bean.isThreadAllocatedMemoryEnabled
    if !enabled then bean.setThreadAllocatedMemoryEnabled(true)
    try
      val fir = Fir.causal(Vec(0.25, 0.5, 0.25)).orThrow
      val runner = fir.newRunner()
      val input = Vec.fill(128)(1.0)
      val output = MutableDVec.zeros(128)
      def calls(count: Int, owned: Boolean): Double =
        var i = 0
        var sum = 0.0
        while i < count do
          if owned then sum += fir.process(input).orThrow(0)
          else
            runner.reset()
            runner.processInto(input, output).orThrow
            sum += output(0)
          i += 1
        sum
      assertEquals(calls(4000, false), 1000.0)
      assertEquals(calls(4000, true), 1000.0)
      val thread = Thread.currentThread().threadId()
      val beforeReuse = bean.getThreadAllocatedBytes(thread)
      val reuseSum = calls(1000, false)
      val reused = bean.getThreadAllocatedBytes(thread) - beforeReuse
      val beforeOwned = bean.getThreadAllocatedBytes(thread)
      val ownedSum = calls(1000, true)
      val owned = bean.getThreadAllocatedBytes(thread) - beforeOwned
      assertEquals(reuseSum, 250.0)
      assertEquals(ownedSum, 250.0)
      println(s"FIR allocation: reused=$reused owned=$owned bytes/1000 calls")
      assert(reused < 128000 && reused < owned / 10,
        s"reused=$reused owned=$owned; a reused runner must not allocate sample/state payloads per call")
    finally if !enabled then bean.setThreadAllocatedMemoryEnabled(false)
