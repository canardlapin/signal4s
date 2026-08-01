package signal4s.bench

import gale.linalg.{DVec, Vec}
import signal4s.*
import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Paths}
import com.sun.management.ThreadMXBean
import java.lang.management.ManagementFactory

/** Lightweight allocation + throughput receipt for the direct convolution operator. */
object DirectOperatorBench:

  def main(args: Array[String]): Unit =
    val kernel = Kernel.causal(Vec.tabulate(64)(i => if i == 0 then 1.0 else 0.0)).orThrow
    val n = 10_000
    val x = DVec.tabulate(n)(i => math.sin(0.01 * i))
    val op = Convolution.operator(kernel, n, OutputRegion.Full).orThrow

    var w = 0
    while w < 20 do
      val _ = op(x)
      w += 1

    val bean = ManagementFactory.getThreadMXBean match
      case sun: ThreadMXBean => Some(sun)
      case _                 => None
    bean.foreach { b =>
      if b.isThreadAllocatedMemorySupported then b.setThreadAllocatedMemoryEnabled(true)
    }
    val threadId = Thread.currentThread().threadId

    val iters = 100
    val t0 = System.nanoTime()
    val alloc0 = bean.map(_.getThreadAllocatedBytes(threadId)).getOrElse(-1L)
    var i = 0
    var sink = 0.0
    while i < iters do
      val y = op(x)
      sink += y(0) + y(y.length - 1)
      i += 1
    val alloc1 = bean.map(_.getThreadAllocatedBytes(threadId)).getOrElse(-1L)
    val t1 = System.nanoTime()

    val ms = (t1 - t0).toDouble / 1e6
    val perCallMs = ms / iters
    val bytesPerCall =
      if alloc0 >= 0 && alloc1 >= alloc0 then (alloc1 - alloc0).toDouble / iters
      else Double.NaN

    val receipt =
      s"""# E3 direct convolution operator bench
         |
         |generated: ${java.time.Instant.now()}
         |inputLength: $n
         |kernelLength: ${kernel.length}
         |region: Full
         |iters: $iters
         |totalMs: ${f"$ms%.3f"}
         |msPerCall: ${f"$perCallMs%.4f"}
         |approxBytesPerCall: ${if bytesPerCall.isNaN then "n/a" else f"$bytesPerCall%.0f"}
         |sink: $sink
         |
         |Notes: allocating `op(x)` path (immutable result). Destination-writing
         |forms arrive with FFT/OLA plans; this receipt baselines the current
         |DirectConvolution-backed operator.
         |""".stripMargin

    val outDir = Paths.get("benchmarks", "receipts")
    Files.createDirectories(outDir)
    val out = outDir.resolve("e3-direct-operator.md")
    Files.writeString(out, receipt, StandardCharsets.UTF_8)
    println(receipt)
    println(s"wrote $out")
