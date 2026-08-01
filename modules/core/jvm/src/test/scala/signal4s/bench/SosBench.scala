package signal4s.bench

import gale.linalg.Vec
import signal4s.*
import signal4s.filter.*
import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Paths}

/** Throughput receipt for SOS filtering (E10 audit). */
object SosBench:
  def main(args: Array[String]): Unit =
    // Butter-like 2-section SOS (from smoke.sosfilt_butter4 shape)
    val sos = SecondOrderCascade
      .fromSosMatrix(
        IArray(
          IArray(0.00482434, 0.00964869, 0.00482434, 1.0, -1.04859958, 0.29614036),
          IArray(1.0, 2.0, 1.0, 1.0, -1.32091343, 0.63273879)
        )
      )
      .orThrow
    val x = Vec.tabulate(8192)(i => math.sin(0.02 * i))
    var w = 0
    while w < 30 do
      val _ = sos.process(x)
      w += 1
    val iters = 100
    val t0 = System.nanoTime()
    var i = 0
    var sink = 0.0
    while i < iters do
      sink += sos.process(x).orThrow(0)
      i += 1
    val ms = (System.nanoTime() - t0).toDouble / 1e6
    val sb = new StringBuilder
    sb.append("# E10 SOS filter bench\n\n")
    sb.append(s"generated: ${java.time.Instant.now()}\n\n")
    sb.append(
      f"n=${x.length} sections=${sos.sectionCount} iters=$iters totalMs=$ms%.3f msPerCall=${ms / iters}%.4f sink=$sink\n"
    )
    val outDir = Paths.get("benchmarks", "receipts")
    Files.createDirectories(outDir)
    val out = outDir.resolve("e10-sos.md")
    Files.writeString(out, sb.toString, StandardCharsets.UTF_8)
    println(sb.toString)
    println(s"wrote $out")
