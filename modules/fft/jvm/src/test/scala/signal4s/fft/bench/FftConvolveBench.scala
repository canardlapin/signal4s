package signal4s.fft.bench

import gale.linalg.Vec
import signal4s.*
import signal4s.fft.FftBackend
import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Paths}

object FftConvolveBench:
  def main(args: Array[String]): Unit =
    FftBackend.ensureInstalled()
    val sb = new StringBuilder
    sb.append("# E6 FFT / OLA convolution bench\n\n")
    sb.append(s"generated: ${java.time.Instant.now()}\n\n")

    val shapes = List((128, 16), (512, 32), (2048, 64), (8192, 128))
    shapes.foreach { case (n, m) =>
      val x = Vec.tabulate(n)(i => math.sin(0.02 * i))
      val k = Kernel.causal(Vec.tabulate(m)(i => 1.0 / (i + 1))).orThrow
      val oneShot = time(40) {
        Convolution(x, k, OutputRegion.Full, ConvolutionMethod.Fft).orThrow(0)
      }
      val plan = Convolution.plan(k, n, OutputRegion.Full, ConvolutionMethod.Fft).orThrow
      val reused = time(40) {
        plan(x).orThrow(0)
      }
      val olaOneShot = time(40) {
        Convolution(x, k, OutputRegion.Full, ConvolutionMethod.OverlapAdd(64)).orThrow(0)
      }
      val olaPlan =
        Convolution.plan(k, n, OutputRegion.Full, ConvolutionMethod.OverlapAdd(64)).orThrow
      val olaReused = time(40) {
        olaPlan(x).orThrow(0)
      }
      sb.append(
        f"N=$n M=$m fftOneShotMs=$oneShot%.4f fftPlanReuseMs=$reused%.4f " +
          f"ola64OneShotMs=$olaOneShot%.4f ola64PlanReuseMs=$olaReused%.4f\n"
      )
    }

    sb.append("\n# Circular convolution\n\n")
    List((64, 8), (128, 16), (256, 16), (512, 32), (1024, 48), (4096, 80), (8192, 128))
      .foreach { case (n, m) =>
      val x = Vec.tabulate(n)(i => math.sin(0.02 * i))
      val k = Kernel.causal(Vec.tabulate(m)(i => 1.0 / (i + 1))).orThrow
      val direct = time(40) {
        Convolution.circular(x, k, n).orThrow(0)
      }
      val fft = time(40) {
        FftBackend.circular(x, k, n).orThrow(0)
      }
      sb.append(f"N=$n M=$m circularDirectMs=$direct%.4f circularFftMs=$fft%.4f\n")
    }

    sb.append("\n# Planned FFT region extraction\n\n")
    List((8192, 128), (16384, 2048)).foreach { case (n, m) =>
      val x = Vec.tabulate(n)(i => math.sin(0.02 * i))
      val k = Kernel.at(Vec.tabulate(m)(i => 1.0 / (i + 1)), zeroLagIndex = m / 2).orThrow
      List(
        "valid" -> OutputRegion.Valid,
        "input" -> OutputRegion.Input(Boundary.Zero)
      ).foreach { case (label, region) =>
        val plan = Convolution.plan(k, n, region, ConvolutionMethod.Fft).orThrow
        val ms = time(40) {
          plan(x).orThrow(0)
        }
        sb.append(f"N=$n M=$m region=$label fftPlanReuseMs=$ms%.4f\n")
      }
    }

    val outDir = Paths.get("benchmarks", "receipts")
    Files.createDirectories(outDir)
    val out = outDir.resolve("e6-fft-convolve.md")
    Files.writeString(out, sb.toString, StandardCharsets.UTF_8)
    println(sb.toString)
    println(s"wrote $out")

  private def time(iters: Int)(body: => Double): Double =
    var w = 0
    while w < 15 do
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
