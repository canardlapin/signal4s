package signal4s.fft.internal

/** JTransforms must remain numerically interchangeable with the portable radix-2 engine. */
class PlatformFftSuite extends munit.FunSuite:

  private val lengths = List(2, 4, 8, 16, 64, 256, 1024, 4096)
  private val realLengths = lengths ++ List(12, 18, 30, 2160)

  test("JVM power-of-two backend matches portable radix-2 forward and inverse"):
    lengths.foreach { n =>
      val re0 = Array.tabulate(n)(i => math.sin(0.17 * i) + 0.03 * (i % 7))
      val im0 = Array.tabulate(n)(i => math.cos(0.11 * i) - 0.02 * (i % 5))

      List(false, true).foreach { inverse =>
        val expectedRe = re0.clone()
        val expectedIm = im0.clone()
        FftEngine.radix2BitrevForTest(expectedRe, expectedIm, inverse)

        val actualRe = re0.clone()
        val actualIm = im0.clone()
        assert(
          PlatformFft.pow2(actualRe, actualIm, inverse),
          clue = s"JTransforms did not accept power-of-two n=$n"
        )

        var i = 0
        while i < n do
          val scale = math.max(
            1.0,
            math.max(
              math.abs(expectedRe(i)),
              math.abs(expectedIm(i))
            )
          )
          assertEqualsDouble(actualRe(i), expectedRe(i), 1e-10 * scale, clue = s"n=$n inv=$inverse re@$i")
          assertEqualsDouble(actualIm(i), expectedIm(i), 1e-10 * scale, clue = s"n=$n inv=$inverse im@$i")
          i += 1
      }
    }

  test("JVM complex backend matches portable mixed-radix forward and inverse"):
    List(6, 12, 60, 2160).foreach { n =>
      val re0 = Array.tabulate(n)(i => math.sin(0.17 * i) + 0.03 * (i % 7))
      val im0 = Array.tabulate(n)(i => math.cos(0.11 * i) - 0.02 * (i % 5))

      List(false, true).foreach { inverse =>
        val expectedRe = re0.clone()
        val expectedIm = im0.clone()
        FftEngine.mixedRadixPortableForTest(expectedRe, expectedIm, inverse)

        val actualRe = re0.clone()
        val actualIm = im0.clone()
        assert(
          PlatformFft.complex(actualRe, actualIm, inverse),
          clue = s"JTransforms did not accept smooth n=$n"
        )

        var i = 0
        while i < n do
          val scale = math.max(
            1.0,
            math.max(
              math.abs(expectedRe(i)),
              math.abs(expectedIm(i))
            )
          )
          assertEqualsDouble(actualRe(i), expectedRe(i), 1e-10 * scale, clue = s"n=$n inv=$inverse re@$i")
          assertEqualsDouble(actualIm(i), expectedIm(i), 1e-10 * scale, clue = s"n=$n inv=$inverse im@$i")
          i += 1
      }
    }

  test("JVM packed-real backend matches portable complex FFT and unscaled inverse"):
    realLengths.foreach { n =>
      val input = Array.tabulate(n)(i => math.sin(0.13 * i) + 0.04 * (i % 5))
      val expectedRe = input.clone()
      val expectedIm = Array.ofDim[Double](n)
      FftEngine.forward(expectedRe, expectedIm)

      val actualRe = input.clone()
      val actualIm = Array.ofDim[Double](n)
      assert(PlatformFft.realForward(actualRe, actualIm, onesided = true))
      var k = 0
      while k <= n / 2 do
        assertEqualsDouble(actualRe(k), expectedRe(k), 1e-10 * n, clue = s"n=$n re@$k")
        assertEqualsDouble(actualIm(k), expectedIm(k), 1e-10 * n, clue = s"n=$n im@$k")
        k += 1

      assert(PlatformFft.realInverse(actualRe, actualIm))
      var i = 0
      while i < n do
        assertEqualsDouble(actualRe(i), n.toDouble * input(i), 1e-10 * n, clue = s"n=$n out@$i")
        assertEqualsDouble(actualIm(i), 0.0, 0.0, clue = s"n=$n imag@$i")
        i += 1
    }
