package signal4s.fft

import gale.linalg.Vec
import signal4s.*
import signal4s.fft.internal.FftEngine

class PackedRealFftSuite extends munit.FunSuite:

  test("forwardReal matches complex forward on real input"):
    List(8, 16, 32, 64, 256).foreach { n =>
      val x = Vec.tabulate(n)(i => math.sin(0.3 * i) + 0.1 * i)
      val reC = Array.ofDim[Double](n)
      val imC = Array.ofDim[Double](n)
      x.copyTo(reC)
      FftEngine.forward(reC, imC)

      val reR = Array.ofDim[Double](n)
      val imR = Array.ofDim[Double](n)
      val packRe = Array.ofDim[Double](n / 2)
      val packIm = Array.ofDim[Double](n / 2)
      x.copyTo(reR)
      FftEngine.forwardRealPortableForTest(reR, imR, packRe, packIm, onesided = false)

      var k = 0
      while k < n do
        assertEqualsDouble(reR(k), reC(k), 1e-9, clue = s"n=$n re@$k")
        assertEqualsDouble(imR(k), imC(k), 1e-9, clue = s"n=$n im@$k")
        k += 1
    }

  test("inverseReal round-trips with forwardReal (unnormalized)"):
    List(8, 16, 32, 128).foreach { n =>
      val x = Vec.tabulate(n)(i => math.cos(0.2 * i) - 0.05 * i)
      val re = Array.ofDim[Double](n)
      val im = Array.ofDim[Double](n)
      val packRe = Array.ofDim[Double](n / 2)
      val packIm = Array.ofDim[Double](n / 2)
      x.copyTo(re)
      FftEngine.forwardRealPortableForTest(re, im, packRe, packIm, onesided = false)
      FftEngine.inverseRealPortableForTest(re, im, packRe, packIm)
      var i = 0
      while i < n do
        assertEqualsDouble(re(i), n.toDouble * x(i), 1e-8, clue = s"n=$n i=$i")
        i += 1
    }

  test("forwardRealOnesided matches full onesided bins; scaled inverse recovers"):
    List(16, 64, 256).foreach { n =>
      val x = Vec.tabulate(n)(i => math.sin(0.13 * i))
      val reFull = Array.ofDim[Double](n)
      val imFull = Array.ofDim[Double](n)
      val reSide = Array.ofDim[Double](n)
      val imSide = Array.ofDim[Double](n)
      val packRe = Array.ofDim[Double](n / 2)
      val packIm = Array.ofDim[Double](n / 2)
      x.copyTo(reFull)
      FftEngine.forwardRealPortableForTest(reFull, imFull, packRe, packIm, onesided = false)
      x.copyTo(reSide)
      FftEngine.forwardRealPortableForTest(reSide, imSide, packRe, packIm, onesided = true)
      val N = n >> 1
      var k = 0
      while k <= N do
        assertEqualsDouble(reSide(k), reFull(k), 1e-9, clue = s"n=$n re@$k")
        assertEqualsDouble(imSide(k), imFull(k), 1e-9, clue = s"n=$n im@$k")
        k += 1
      val out = Array.ofDim[Double](n)
      FftEngine.inverseRealScaledToPortableForTest(
        reSide,
        imSide,
        packRe,
        packIm,
        1.0 / n,
        out,
        n
      )
      k = 0
      while k < n do
        assertEqualsDouble(out(k), x(k), 1e-8, clue = s"n=$n out@$k")
        k += 1
    }

  test("FFT convolution plan matches direct"):
    FftBackend.ensureInstalled()
    val n = 64
    val m = 9
    val x = Vec.tabulate(n)(i => math.sin(0.15 * i))
    val k = Kernel.causal(Vec.tabulate(m)(i => 1.0 / (i + 1.0))).orThrow
    val direct = Convolution(x, k, OutputRegion.Full, ConvolutionMethod.Direct).orThrow
    val plan = Convolution.plan(k, n, OutputRegion.Full, ConvolutionMethod.Fft).orThrow
    val fft = plan(x).orThrow
    assertEquals(fft.length, direct.length)
    var i = 0
    while i < fft.length do
      assertEqualsDouble(fft(i), direct(i), 1e-9, clue = s"i=$i")
      i += 1
