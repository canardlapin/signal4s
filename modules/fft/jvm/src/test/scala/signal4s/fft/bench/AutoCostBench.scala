package signal4s.fft.bench

import gale.linalg.Vec
import signal4s.*
import signal4s.fft.FftBackend
import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Paths}

/** Microbench used to calibrate [[signal4s.fft.AutoCostModel]] coefficients. */
object AutoCostBench:
  def main(args: Array[String]): Unit =
    FftBackend.ensureInstalled()
    val cases = List(
      (64, 8),
      (256, 16),
      (1024, 32),
      (4096, 64),
      (8192, 9),
      (1024, 512),
      (8192, 128),
      (65536, 64),
      (16384, 2048),
      (16384, 4096)
    )
    val sb = new StringBuilder
    sb.append("# E6 Auto cost model calibration\n\n")
    sb.append(s"generated: ${java.time.Instant.now()}\n\n")
    sb.append("Planned wall-clock (ms/call) for Full convolution; warm runs.\n\n")
    sb.append("| N | M | Direct | FFT | OLA(B) | Auto picks |\n")
    sb.append("|---|---|--------|-----|--------|------------|\n")

    cases.foreach { case (n, m) =>
      val x = Vec.tabulate(n)(i => math.sin(0.01 * i))
      val k = Kernel.causal(Vec.tabulate(m)(i => 0.1 * (i + 1.0))).orThrow
      val block = signal4s.fft.AutoCostModel.defaultBlock(n, m)
      val directPlan =
        Convolution.plan(k, n, OutputRegion.Full, ConvolutionMethod.Direct).orThrow
      val fftPlan =
        Convolution.plan(k, n, OutputRegion.Full, ConvolutionMethod.Fft).orThrow
      val olaPlan =
        Convolution.plan(k, n, OutputRegion.Full, ConvolutionMethod.OverlapAdd(block)).orThrow
      val d = timeMs(80) {
        directPlan(x).orThrow(0)
      }
      val f = timeMs(80) {
        fftPlan(x).orThrow(0)
      }
      val o = timeMs(80) {
        olaPlan(x).orThrow(0)
      }
      val plan = Convolution.plan(k, n, OutputRegion.Full, ConvolutionMethod.Auto).orThrow
      sb.append(
        f"| $n | $m | $d%.4f | $f%.4f | $o%.4f (B=$block) | ${plan.selectedMethod} |\n"
      )
    }
    sb.append("\n## Model\n\n")
    sb.append("AutoCostModel uses relative costs:\n")
    sb.append("- Direct ≈ 1.0 · N · M\n")
    sb.append("- FFT ≈ 12.0 · L · log2(L), L = nextPow2(N+M-1)\n")
    sb.append("- OLA ≈ 12.0 · ⌈N/B⌉ · L_b · log2(L_b), L_b = nextPow2(B+M-1)\n")
    sb.append("\nCoefficients chosen so ordering tracks this JVM portable matrix.\n")

    val outDir = Paths.get("benchmarks", "receipts")
    Files.createDirectories(outDir)
    val out = outDir.resolve("e6-auto-cost.md")
    Files.writeString(out, sb.toString, StandardCharsets.UTF_8)
    println(sb.toString)
    println(s"wrote $out")

  private def timeMs(iters: Int)(body: => Double): Double =
    var w = 0
    while w < 20 do
      val _ = body
      w += 1
    val t0 = System.nanoTime()
    var i = 0
    var sink = 0.0
    while i < iters do
      sink += body
      i += 1
    val _ = sink
    (System.nanoTime() - t0).toDouble / 1e6 / iters
