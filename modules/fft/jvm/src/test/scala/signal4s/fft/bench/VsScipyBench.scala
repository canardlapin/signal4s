package signal4s.fft.bench

import gale.linalg.{DVec, MutableDVec, Vec}
import signal4s.*
import signal4s.filter.SecondOrderCascade
import signal4s.fft.{ComplexVector, FftBackend, FftNormalization, FftPlan}
import signal4s.multirate.Upfirdn
import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Paths}

/** Paired JVM timings for docs/PERFORMANCE.md vs-SciPy receipt.
  *
  * Shape grids match `fixtures/generate/bench_vs_scipy.py` — broad coverage, not
  * a shortlist chosen to flatter a single optimization.
  */
object VsScipyBench:

  private val Warmup = 50
  private val Iters = 100
  private val Trials = 7

  def main(args: Array[String]): Unit =
    FftBackend.ensureInstalled()
    // First full pass settles JIT/GC; receipt uses the second pass only.
    val _ = measureAll()
    val rows = measureAll()

    val scalaJson =
      ujson(
        "generated" -> java.time.Instant.now().toString,
        "warmup" -> Warmup,
        "iters" -> Iters,
        "trials" -> Trials,
        "rows" -> rows
      )
    val outDir = Paths.get("benchmarks", "receipts")
    Files.createDirectories(outDir)
    val scalaPath = outDir.resolve("vs-scipy-scala.json")
    Files.writeString(scalaPath, scalaJson, StandardCharsets.UTF_8)
    println(s"wrote $scalaPath")

    val pyPath = outDir.resolve("vs-scipy-python.json")
    val md = mergeMarkdown(pyPath, rows)
    val mdPath = outDir.resolve("vs-scipy.md")
    Files.writeString(mdPath, md, StandardCharsets.UTF_8)
    println(md)
    println(s"wrote $mdPath")

  private def measureAll(): List[Row] =
    convolveShapes.flatMap { case (n, m) =>
      List(convolveDirect(n, m), convolveFftOneShot(n, m), convolveFftPlanned(n, m))
    } ++
      causalFirShapes.flatMap { case (n, m) =>
        List(causalFir(n, m), convolveValidDirect(n, m), circularFft(n, m))
      } ++
      fftLengths.flatMap(n => List(fftPow2(n), rfftPow2(n))) ++
      smoothFftLengths.map(fftSmooth) ++
      sosLengths.map(n => sosfilt(n, 2)) ++
      upfirdnShapes.map { case (n, m, up, down) => upfirdn(n, m, up, down) }

  private val convolveShapes: List[(Int, Int)] =
    List(
      (128, 9),
      (256, 16),
      (512, 33),
      (1024, 48),
      (2048, 64),
      (4096, 80),
      (8192, 128),
      (3000, 50)
    )

  private val fftLengths: List[Int] =
    List(64, 256, 1024, 2048, 4096, 8192)

  private val smoothFftLengths: List[Int] =
    List(12, 60, 2160, 8640)

  private val causalFirShapes: List[(Int, Int)] =
    List((1024, 48), (4096, 80), (8192, 128))

  private val sosLengths: List[Int] =
    List(512, 2048, 8192, 32768, 65536)

  private val upfirdnShapes: List[(Int, Int, Int, Int)] =
    List(
      (512, 15, 2, 3),
      (2048, 31, 3, 2),
      (4096, 47, 4, 3),
      (8192, 63, 5, 3),
      (16384, 63, 2, 1),
      (8192, 63, 1, 2)
    )

  private def convolveDirect(n: Int, m: Int): Row =
    val x = Vec.tabulate(n)(i => math.sin(0.02 * i))
    val k = Kernel.causal(Vec.tabulate(m)(i => 1.0 / (i + 1))).orThrow
    val plan = Convolution.plan(k, n, OutputRegion.Full, ConvolutionMethod.Direct).orThrow
    val ms = timeMs {
      plan(x).orThrow(0)
    }
    Row(s"convolve_direct_full_n${n}_m$m", "convolve_direct_full", ms, Map("n" -> n, "m" -> m))

  private def convolveFftOneShot(n: Int, m: Int): Row =
    val x = Vec.tabulate(n)(i => math.sin(0.02 * i))
    val k = Kernel.causal(Vec.tabulate(m)(i => 1.0 / (i + 1))).orThrow
    val ms = timeMs {
      Convolution(x, k, OutputRegion.Full, ConvolutionMethod.Fft).orThrow(0)
    }
    Row(s"convolve_fft_full_n${n}_m$m", "convolve_fft_full", ms, Map("n" -> n, "m" -> m))

  private def convolveFftPlanned(n: Int, m: Int): Row =
    val x = Vec.tabulate(n)(i => math.sin(0.02 * i))
    val k = Kernel.causal(Vec.tabulate(m)(i => 1.0 / (i + 1))).orThrow
    val plan = Convolution.plan(k, n, OutputRegion.Full, ConvolutionMethod.Fft).orThrow
    val ms = timeMs {
      plan(x).orThrow(0)
    }
    Row(
      s"convolve_fft_planned_full_n${n}_m$m",
      "convolve_fft_planned_full",
      ms,
      Map("n" -> n, "m" -> m)
    )

  private def causalFir(n: Int, m: Int): Row =
    val x = Vec.tabulate(n)(i => math.sin(0.02 * i))
    val k = Kernel.causal(Vec.tabulate(m)(i => 1.0 / (i + 1))).orThrow
    val plan =
      Convolution.plan(k, n, OutputRegion.Input(Boundary.Zero), ConvolutionMethod.Direct).orThrow
    val ms = timeMs {
      plan(x).orThrow(0)
    }
    Row(s"fir_causal_n${n}_m$m", "fir_causal", ms, Map("n" -> n, "m" -> m))

  private def circularFft(n: Int, m: Int): Row =
    val x = Vec.tabulate(n)(i => math.sin(0.02 * i))
    val k = Kernel.causal(Vec.tabulate(m)(i => 1.0 / (i + 1))).orThrow
    val ms = timeMs {
      Convolution.circular(x, k, n).orThrow(0)
    }
    Row(s"circular_fft_n${n}_m$m", "circular_fft", ms, Map("n" -> n, "m" -> m))

  private def convolveValidDirect(n: Int, m: Int): Row =
    val x = Vec.tabulate(n)(i => math.sin(0.02 * i))
    val k = Kernel.causal(Vec.tabulate(m)(i => 1.0 / (i + 1))).orThrow
    val plan = Convolution.plan(k, n, OutputRegion.Valid, ConvolutionMethod.Direct).orThrow
    val ms = timeMs {
      plan(x).orThrow(0)
    }
    Row(
      s"convolve_valid_direct_n${n}_m$m",
      "convolve_valid_direct",
      ms,
      Map("n" -> n, "m" -> m)
    )

  private def fftPow2(n: Int): Row =
    val re = Vec.tabulate(n)(i => math.sin(0.02 * i))
    val im = DVec.zeros(n)
    val cv = ComplexVector(re, im).orThrow
    val plan = FftPlan(n, FftNormalization.Backward).orThrow
    val ws = plan.newWorkspace()
    val ms = timeMs {
      plan.forwardInto(cv, ws).orThrow
      ws.length.toDouble
    }
    Row(s"fft_pow2_n$n", "fft_pow2", ms, Map("n" -> n))

  private def rfftPow2(n: Int): Row =
    val x = Vec.tabulate(n)(i => math.sin(0.02 * i))
    val plan = signal4s.fft.RealFftPlan(n, FftNormalization.Backward).orThrow
    val ws = plan.newWorkspace()
    val ms = timeMs {
      plan.forwardInto(x, ws).orThrow
      ws.length.toDouble
    }
    Row(s"rfft_pow2_n$n", "rfft_pow2", ms, Map("n" -> n))

  private def fftSmooth(n: Int): Row =
    val re = Vec.tabulate(n)(i => math.sin(0.02 * i))
    val im = DVec.tabulate(n)(i => math.cos(0.017 * i))
    val cv = ComplexVector(re, im).orThrow
    val plan = FftPlan(n, FftNormalization.Backward).orThrow
    val ws = plan.newWorkspace()
    val ms = timeMs {
      plan.forwardInto(cv, ws).orThrow
      ws.length.toDouble
    }
    Row(s"fft_smooth_n$n", "fft_smooth", ms, Map("n" -> n))

  private def sosfilt(n: Int, sections: Int): Row =
    val sos = SecondOrderCascade
      .fromSosMatrix(
        IArray(
          IArray(0.00482434, 0.00964869, 0.00482434, 1.0, -1.04859958, 0.29614036),
          IArray(1.0, 2.0, 1.0, 1.0, -1.32091343, 0.63273879)
        )
      )
      .orThrow
    require(sos.sectionCount == sections)
    val x = Vec.tabulate(n)(i => math.sin(0.02 * i))
    val runner = sos.newRunner()
    val out = MutableDVec.zeros(n)
    val ms = timeMs {
      runner.reset()
      runner.processInto(x, out).orThrow
      out(0)
    }
    Row(s"sosfilt_n${n}_sec$sections", "sosfilt", ms, Map("n" -> n, "sections" -> sections))

  private def upfirdn(n: Int, m: Int, up: Int, down: Int): Row =
    val x = Vec.tabulate(n)(i => math.sin(0.02 * i))
    val h = Vec.tabulate(m) { i =>
      val t = if m == 1 then 0.5 else i.toDouble / (m - 1)
      0.54 - 0.46 * math.cos(2.0 * math.Pi * t)
    }
    val ms = timeMs {
      Upfirdn(h, x, up, down).orThrow(0)
    }
    Row(
      s"upfirdn_n${n}_m${m}_up${up}_down$down",
      "upfirdn",
      ms,
      Map("n" -> n, "m" -> m, "up" -> up, "down" -> down)
    )

  private def timeMs(body: => Double): Double =
    var w = 0
    while w < Warmup do
      val _ = body
      w += 1
    val samples = Array.ofDim[Double](Trials)
    var t = 0
    while t < Trials do
      val t0 = System.nanoTime()
      var i = 0
      var sink = 0.0
      while i < Iters do
        sink += body
        i += 1
      val _ = sink
      samples(t) = (System.nanoTime() - t0).toDouble / 1e6 / Iters
      t += 1
    java.util.Arrays.sort(samples)
    samples(Trials / 2)

  private def mergeMarkdown(pyPath: java.nio.file.Path, scalaRows: List[Row]): String =
    val sb = new StringBuilder
    sb.append("# Performance vs SciPy\n\n")
    sb.append(s"generated: ${java.time.Instant.now()}\n\n")
    sb.append("See [`docs/PERFORMANCE.md`](../../docs/PERFORMANCE.md) for fairness rules.\n\n")
    sb.append(s"Trials: median of $Trials × $Iters iters after $Warmup warm-up.\n\n")
    if !Files.exists(pyPath) then
      sb.append(
        s"_Python receipt missing at `$pyPath`. Run `python3 fixtures/generate/bench_vs_scipy.py` first._\n"
      )
      scalaRows.foreach(r => sb.append(f"- `${r.id}`: scala=${r.ms}%.4f ms\n"))
      return sb.toString

    val pyText = Files.readString(pyPath)
    val pyById = parsePythonMs(pyText)
    sb.append("| id | scipy ms | signal4s ms | ratio (s4s/scipy) |\n")
    sb.append("|---|---:|---:|---:|\n")
    scalaRows.foreach { r =>
      pyById.get(r.id) match
        case Some(pyMs) =>
          val ratio = r.ms / pyMs
          sb.append(f"| `${r.id}` | $pyMs%.4f | ${r.ms}%.4f | $ratio%.2f× |\n")
        case None =>
          sb.append(f"| `${r.id}` | — | ${r.ms}%.4f | — |\n")
    }
    sb.toString

  private def parsePythonMs(text: String): Map[String, Double] =
    val rowPat =
      """\{[^{}]*"id"\s*:\s*"([^"]+)"[^{}]*"ms"\s*:\s*([0-9.eE+-]+)[^{}]*\}|\{[^{}]*"ms"\s*:\s*([0-9.eE+-]+)[^{}]*"id"\s*:\s*"([^"]+)"[^{}]*\}""".r
    rowPat
      .findAllMatchIn(text)
      .map { m =>
        if m.group(1) != null then m.group(1) -> m.group(2).toDouble
        else m.group(4) -> m.group(3).toDouble
      }
      .toMap

  private final case class Row(id: String, op: String, ms: Double, dims: Map[String, Int]):
    def json: String =
      val dimFields = dims.toList
        .sortBy(_._1)
        .map { case (k, v) => s""""$k": $v""" }
        .mkString(", ")
      s"""{"id":"$id","op":"$op","ms":$ms${if dimFields.isEmpty then "" else ", " + dimFields}}"""

  private def ujson(fields: (String, Any)*): String =
    val body = fields
      .map {
        case (k, v: String)                   => s""""$k":"$v""""
        case (k, v: Int)                      => s""""$k":$v"""
        case (k, v: List[Row] @unchecked)     => s""""$k":[${v.map(_.json).mkString(",")}]"""
        case (k, v)                           => s""""$k":"$v""""
      }
      .mkString(",\n  ")
    s"{\n  $body\n}\n"
