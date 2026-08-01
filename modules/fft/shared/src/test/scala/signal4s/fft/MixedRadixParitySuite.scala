package signal4s.fft

import signal4s.fft.internal.FftEngine

/** Mixed-radix lengths must match the unscaled DFT (not just pow2 / fixture cells). */
class MixedRadixParitySuite extends munit.FunSuite:

  private val lengths =
    List(6, 12, 15, 18, 20, 24, 30, 36, 45, 48, 60, 72, 90, 120, 180, 240, 360, 720, 1080, 2160)

  test("forward mixed-radix matches naive DFT"):
    lengths.foreach { n =>
      val re = Array.tabulate(n)(i => math.sin(0.31 * i) + 0.02 * (i % 7))
      val im = Array.tabulate(n)(i => math.cos(0.17 * i) - 0.01 * (i % 5))
      val expRe = Array.ofDim[Double](n)
      val expIm = Array.ofDim[Double](n)
      naiveDft(re, im, expRe, expIm, inverse = false)
      val gotRe = re.clone()
      val gotIm = im.clone()
      FftEngine.mixedRadixPortableForTest(gotRe, gotIm, inverse = false)
      assertClose(gotRe, gotIm, expRe, expIm, clue = s"forward n=$n")
    }

  test("inverse mixed-radix round-trips (unscaled)"):
    lengths.foreach { n =>
      val re = Array.tabulate(n)(i => 0.1 * i - 0.03 * (i % 11))
      val im = Array.tabulate(n)(i => math.sin(0.05 * i))
      val origRe = re.clone()
      val origIm = im.clone()
      FftEngine.mixedRadixPortableForTest(re, im, inverse = false)
      FftEngine.mixedRadixPortableForTest(re, im, inverse = true)
      var i = 0
      while i < n do
        val tol = 1e-9 + 1e-9 * n * math.max(math.abs(origRe(i)), math.abs(origIm(i)))
        assert(math.abs(re(i) - origRe(i) * n) <= tol, clue = s"n=$n i=$i re")
        assert(math.abs(im(i) - origIm(i) * n) <= tol, clue = s"n=$n i=$i im")
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

  private def assertClose(
      gotRe: Array[Double],
      gotIm: Array[Double],
      expRe: Array[Double],
      expIm: Array[Double],
      clue: String
  ): Unit =
    var i = 0
    while i < gotRe.length do
      val scale = math.max(math.abs(expRe(i)), math.abs(expIm(i)))
      val tol = 1e-9 + 1e-9 * gotRe.length * scale
      val err = math.hypot(gotRe(i) - expRe(i), gotIm(i) - expIm(i))
      assert(err <= tol, clue = s"$clue@$i err=$err tol=$tol")
      i += 1
