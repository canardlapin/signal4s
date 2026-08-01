package signal4s.fft.bench

import gale.linalg.{DVecBuilder, Vec}
import signal4s.*
import signal4s.fft.internal.{FastFftLength, FftEngine}

/** A1/B1 breakdown: Direct scatter vs adopt; FFT-conv forward/mul/inv/adopt.
  *
  * Not a receipt cell — engine-wide stage costs to guide Phase A/B work.
  */
object HotPathProfile:

  private val Warmup = 60
  private val Iters = 120
  private val Trials = 7

  def main(args: Array[String]): Unit =
    println("=== A1 Direct: copy / scatter / adopt ===")
    List((2048, 64), (4096, 80), (8192, 128), (3000, 50)).foreach { case (n, m) =>
      profileDirect(n, m)
    }
    println()
    println("=== B1 FFT-conv: pad / fwd / mul / inv / adopt ===")
    List((2048, 64), (4096, 80), (8192, 128), (3000, 50)).foreach { case (n, m) =>
      profileFftConv(n, m)
    }

  private def profileDirect(n: Int, m: Int): Unit =
    val x = Vec.tabulate(n)(i => math.sin(0.02 * i))
    val h = Vec.tabulate(m)(i => 1.0 / (i + 1))
    val xa = Array.ofDim[Double](n)
    val ha = Array.ofDim[Double](m)
    val outLen = n + m - 1
    val out = Array.ofDim[Double](outLen)

    def copyIn(): Double =
      x.copyTo(xa)
      h.copyTo(ha)
      xa(0) + ha(0)

    def scatter(): Double =
      java.util.Arrays.fill(out, 0.0)
      scatterFullUnrolled(xa, ha, out)
      out(0)

    def adopt(): Double =
      val b = DVecBuilder.zeros(outLen)
      var i = 0
      while i < outLen do
        b(i) = out(i)
        i += 1
      b.result()(0)

    def endToEnd(): Double =
      val _ = copyIn()
      val _ = scatter()
      adopt()

    val tCopy = medianMs(copyIn())
    val tScatter = medianMs(scatter())
    val tAdopt = medianMs(adopt())
    val tAll = medianMs(endToEnd())
    println(
      f"n=$n%5d m=$m%4d  copy=$tCopy%.4f  scatter=$tScatter%.4f  adopt=$tAdopt%.4f  e2e=$tAll%.4f ms  " +
        f"adoptPct=${100 * tAdopt / tAll}%.0f scatterPct=${100 * tScatter / tAll}%.0f"
    )

  private def profileFftConv(n: Int, m: Int): Unit =
    val x = Vec.tabulate(n)(i => math.sin(0.02 * i))
    val k = Kernel.causal(Vec.tabulate(m)(i => 1.0 / (i + 1))).orThrow
    val outLen = n + m - 1
    val nfft = FastFftLength(outLen)
    val workRe = Array.ofDim[Double](nfft)
    val workIm = Array.ofDim[Double](nfft)
    val packRe = Array.ofDim[Double](nfft / 2)
    val packIm = Array.ofDim[Double](nfft / 2)
    val outBuf = Array.ofDim[Double](outLen)
    val halfRe = Array.ofDim[Double](nfft / 2 + 1)
    val halfIm = Array.ofDim[Double](nfft / 2 + 1)

    k.taps.copyTo(workRe)
    java.util.Arrays.fill(workRe, m, nfft, 0.0)
    FftEngine.forwardRealOnesided(workRe, workIm, packRe, packIm)
    var b = 0
    while b <= nfft / 2 do
      halfRe(b) = workRe(b)
      halfIm(b) = workIm(b)
      b += 1

    def padCopy(): Double =
      x.copyTo(workRe)
      if n < nfft then java.util.Arrays.fill(workRe, n, nfft, 0.0)
      workRe(0)

    def isoFwd(): Double =
      x.copyTo(workRe)
      if n < nfft then java.util.Arrays.fill(workRe, n, nfft, 0.0)
      val t0 = System.nanoTime()
      FftEngine.forwardRealOnesided(workRe, workIm, packRe, packIm)
      (System.nanoTime() - t0).toDouble

    def isoMul(): Double =
      x.copyTo(workRe)
      if n < nfft then java.util.Arrays.fill(workRe, n, nfft, 0.0)
      FftEngine.forwardRealOnesided(workRe, workIm, packRe, packIm)
      val t0 = System.nanoTime()
      mulOnesided(workRe, workIm, halfRe, halfIm, nfft)
      (System.nanoTime() - t0).toDouble

    def isoInv(): Double =
      x.copyTo(workRe)
      if n < nfft then java.util.Arrays.fill(workRe, n, nfft, 0.0)
      FftEngine.forwardRealOnesided(workRe, workIm, packRe, packIm)
      mulOnesided(workRe, workIm, halfRe, halfIm, nfft)
      val t0 = System.nanoTime()
      FftEngine.inverseRealScaledTo(workRe, workIm, packRe, packIm, 1.0 / nfft, outBuf, outLen)
      (System.nanoTime() - t0).toDouble

    def adopt(): Double =
      val builder = DVecBuilder.zeros(outLen)
      var i = 0
      while i < outLen do
        builder(i) = outBuf(i)
        i += 1
      builder.result()(0)

    def endToEnd(): Double =
      x.copyTo(workRe)
      if n < nfft then java.util.Arrays.fill(workRe, n, nfft, 0.0)
      FftEngine.forwardRealOnesided(workRe, workIm, packRe, packIm)
      mulOnesided(workRe, workIm, halfRe, halfIm, nfft)
      FftEngine.inverseRealScaledTo(workRe, workIm, packRe, packIm, 1.0 / nfft, outBuf, outLen)
      adopt()

    val tPad = medianMs(padCopy())
    val tFwd = medianNs(isoFwd()) / 1e6
    val tMul = medianNs(isoMul()) / 1e6
    val tInv = medianNs(isoInv()) / 1e6
    val tAdopt = medianMs(adopt())
    val tE2e = medianMs(endToEnd())
    val stages = tPad + tFwd + tMul + tInv + tAdopt
    println(
      f"n=$n%5d m=$m%4d nfft=$nfft%5d  pad=$tPad%.4f  fwd=$tFwd%.4f  mul=$tMul%.4f  inv=$tInv%.4f  adopt=$tAdopt%.4f  sum=$stages%.4f  e2e=$tE2e%.4f ms  " +
        f"fftPct=${100 * (tFwd + tInv) / stages}%.0f adoptPct=${100 * tAdopt / stages}%.0f"
    )

  /** Mirror of JVM PlatformDirect.scatterFullUnrolled for package-local profiling. */
  private def scatterFullUnrolled(x: Array[Double], h: Array[Double], out: Array[Double]): Unit =
    val n = x.length
    val m = h.length
    val m8 = m & ~7
    var i = 0
    while i < n do
      val xi = x(i)
      if xi != 0.0 then
        var j = 0
        while j < m8 do
          val base = i + j
          out(base) += xi * h(j)
          out(base + 1) += xi * h(j + 1)
          out(base + 2) += xi * h(j + 2)
          out(base + 3) += xi * h(j + 3)
          out(base + 4) += xi * h(j + 4)
          out(base + 5) += xi * h(j + 5)
          out(base + 6) += xi * h(j + 6)
          out(base + 7) += xi * h(j + 7)
          j += 8
        while j < m do
          out(i + j) += xi * h(j)
          j += 1
      i += 1

  private def mulOnesided(
      workRe: Array[Double],
      workIm: Array[Double],
      halfRe: Array[Double],
      halfIm: Array[Double],
      nfft: Int
  ): Unit =
    val N = nfft >> 1
    workRe(0) *= halfRe(0)
    workIm(0) = 0.0
    workRe(N) *= halfRe(N)
    workIm(N) = 0.0
    var kk = 1
    while kk < N do
      val ar = workRe(kk)
      val ai = workIm(kk)
      workRe(kk) = ar * halfRe(kk) - ai * halfIm(kk)
      workIm(kk) = ar * halfIm(kk) + ai * halfRe(kk)
      kk += 1

  private def medianMs(body: => Double): Double =
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

  private def medianNs(body: => Double): Double =
    var w = 0
    while w < Warmup do
      val _ = body
      w += 1
    val samples = Array.ofDim[Double](Trials)
    var t = 0
    while t < Trials do
      var s = 0.0
      var i = 0
      while i < Iters do
        s += body
        i += 1
      samples(t) = s / Iters
      t += 1
    java.util.Arrays.sort(samples)
    samples(Trials / 2)
