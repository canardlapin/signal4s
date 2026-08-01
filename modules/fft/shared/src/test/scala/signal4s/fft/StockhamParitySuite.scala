package signal4s.fft

import signal4s.fft.internal.FftEngine

/** Stockham must match bit-rev radix-2 and the unscaled DFT. */
class StockhamParitySuite extends munit.FunSuite:

  private val lengths = List(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096)

  test("Stockham matches bit-rev radix-2 (forward and inverse)"):
    lengths.foreach { n =>
      val re0 = Array.tabulate(n)(i => math.sin(0.21 * i) + 0.03 * (i % 9))
      val im0 = Array.tabulate(n)(i => math.cos(0.13 * i) - 0.02 * (i % 7))
      List(false, true).foreach { inverse =>
        val bitRe = re0.clone()
        val bitIm = im0.clone()
        FftEngine.radix2BitrevForTest(bitRe, bitIm, inverse)
        val stRe = re0.clone()
        val stIm = im0.clone()
        val tmpRe = Array.ofDim[Double](n)
        val tmpIm = Array.ofDim[Double](n)
        FftEngine.radix2StockhamForTest(stRe, stIm, tmpRe, tmpIm, inverse)
        var i = 0
        while i < n do
          val tol = 1e-9 + 1e-9 * n * math.max(math.abs(bitRe(i)), math.abs(bitIm(i)))
          assert(math.abs(stRe(i) - bitRe(i)) <= tol, clue = s"n=$n inv=$inverse re@$i")
          assert(math.abs(stIm(i) - bitIm(i)) <= tol, clue = s"n=$n inv=$inverse im@$i")
          i += 1
      }
    }

  test("Stockham forward matches naive DFT on medium sizes"):
    List(8, 32, 128, 256).foreach { n =>
      val re = Array.tabulate(n)(i => math.sin(0.31 * i))
      val im = Array.tabulate(n)(i => math.cos(0.17 * i))
      val expRe = Array.ofDim[Double](n)
      val expIm = Array.ofDim[Double](n)
      naiveDft(re, im, expRe, expIm, inverse = false)
      val gotRe = re.clone()
      val gotIm = im.clone()
      val tmpRe = Array.ofDim[Double](n)
      val tmpIm = Array.ofDim[Double](n)
      FftEngine.radix2StockhamForTest(gotRe, gotIm, tmpRe, tmpIm, inverse = false)
      var i = 0
      while i < n do
        val tol = 1e-9 + 1e-9 * n * math.max(math.abs(expRe(i)), math.abs(expIm(i)))
        val err = math.hypot(gotRe(i) - expRe(i), gotIm(i) - expIm(i))
        assert(err <= tol, clue = s"n=$n@$i err=$err")
        i += 1
    }

  private def naiveDft(
      re: Array[Double],
      im: Array[Double],
      outRe: Array[Double],
      outIm: Array[Double],
      inverse: Boolean
  ): Unit =
    val n = re.length
    val sign = if inverse then 1.0 else -1.0
    var k = 0
    while k < n do
      var sr = 0.0
      var si = 0.0
      var t = 0
      while t < n do
        val a = sign * 2.0 * math.Pi * k.toDouble * t.toDouble / n.toDouble
        val wr = math.cos(a)
        val wi = math.sin(a)
        sr += re(t) * wr - im(t) * wi
        si += re(t) * wi + im(t) * wr
        t += 1
      outRe(k) = sr
      outIm(k) = si
      k += 1
