package signal4s.fft.internal

/** Scala.js intentionally uses the portable radix-2 implementation for power-of-two FFTs. */
class PlatformFftSuite extends munit.FunSuite:

  test("fallback declines acceleration and production engine matches radix-2 reference"):
    List(2, 8, 64, 256).foreach { n =>
      val re0 = Array.tabulate(n)(i => math.sin(0.13 * i))
      val im0 = Array.tabulate(n)(i => math.cos(0.07 * i))
      assert(!PlatformFft.pow2(re0.clone(), im0.clone(), inverse = false), clue = s"n=$n")

      List(false, true).foreach { inverse =>
        val expectedRe = re0.clone()
        val expectedIm = im0.clone()
        FftEngine.radix2BitrevForTest(expectedRe, expectedIm, inverse)

        val actualRe = re0.clone()
        val actualIm = im0.clone()
        if inverse then FftEngine.inverse(actualRe, actualIm)
        else FftEngine.forward(actualRe, actualIm)

        var i = 0
        while i < n do
          assertEqualsDouble(actualRe(i), expectedRe(i), 1e-11, clue = s"n=$n inv=$inverse re@$i")
          assertEqualsDouble(actualIm(i), expectedIm(i), 1e-11, clue = s"n=$n inv=$inverse im@$i")
          i += 1
      }
    }
