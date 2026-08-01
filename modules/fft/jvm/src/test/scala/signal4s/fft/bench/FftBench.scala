package signal4s.fft.bench

import gale.linalg.Vec
import signal4s.*
import signal4s.fft.*
import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Paths}

object FftBench:
  def main(args: Array[String]): Unit =
    val lengths = List(64, 256, 1024, 2160, 4096, 8640, 16384, 17, 12)
    val sb = new StringBuilder
    sb.append("# E5 FFT bench\n\n")
    sb.append(s"generated: ${java.time.Instant.now()}\n\n")
    lengths.foreach { n =>
      val plan = RealFftPlan(n, FftNormalization.Backward).orThrow
      val x = Vec.tabulate(n)(i => math.sin(0.01 * i))
      var w = 0
      while w < 50 do
        val _ = plan.forward(x)
        w += 1
      val iters = 200
      val t0 = System.nanoTime()
      var i = 0
      var sink = 0.0
      while i < iters do
        sink += plan.forward(x).orThrow.bins.real(0)
        i += 1
      val ms = (System.nanoTime() - t0).toDouble / 1e6
      sb.append(s"length=$n iters=$iters totalMs=${f"$ms%.3f"} msPerCall=${f"${ms / iters}%.4f"} sink=$sink\n")
    }
    sb.append("\n# Complex non-power-of-two workspace FFT\n\n")
    List(12, 60, 2160, 8640).foreach { n =>
      val plan = FftPlan(n, FftNormalization.Backward).orThrow
      val input = ComplexVector.tabulate(n)(i =>
        Complex(math.sin(0.01 * i), math.cos(0.017 * i))
      ).orThrow
      val ws = plan.newWorkspace()
      var w = 0
      while w < 50 do
        plan.forwardInto(input, ws).orThrow
        w += 1
      val iters = 200
      val t0 = System.nanoTime()
      var i = 0
      var sink = 0.0
      while i < iters do
        plan.forwardInto(input, ws).orThrow
        sink += ws.re(0)
        i += 1
      val ms = (System.nanoTime() - t0).toDouble / 1e6
      sb.append(
        s"complexLength=$n iters=$iters totalMs=${f"$ms%.3f"} msPerCall=${f"${ms / iters}%.4f"} sink=$sink\n"
      )
    }
    val outDir = Paths.get("benchmarks", "receipts")
    Files.createDirectories(outDir)
    val out = outDir.resolve("e5-fft.md")
    Files.writeString(out, sb.toString, StandardCharsets.UTF_8)
    println(sb.toString)
    println(s"wrote $out")
