package signal4s.fft

import signal4s.fft.internal.FastFftLength

class FastFftLengthSuite extends munit.FunSuite:

  test("next power of two"):
    assertEquals(FastFftLength.nextPowerOfTwo(1), 1)
    assertEquals(FastFftLength.nextPowerOfTwo(2), 2)
    assertEquals(FastFftLength.nextPowerOfTwo(3), 4)
    assertEquals(FastFftLength.nextPowerOfTwo(2111), 4096)

  test("next five-smooth is even and 5-smooth"):
    List(1, 7, 100, 125, 2111, 8319, 5000).foreach { min =>
      val n = FastFftLength.nextFiveSmooth(min)
      assert(n >= min, clue = s"min=$min n=$n")
      assert(n % 2 == 0, clue = s"n=$n")
      assert(FastFftLength.isFiveSmooth(n), clue = s"n=$n")
    }

  test("fftLength is power of two ≥ min (packed real path)"):
    List(8, 100, 271, 2111, 4096, 8319).foreach { min =>
      val n = FastFftLength(min)
      assert(n >= min)
      assert(FastFftLength.isPowerOfTwo(n))
    }

  test("five-smooth helpers stay available for non-pow2 engine paths"):
    assertEquals(FastFftLength.nextFiveSmooth(2111), 2160)
    assertEquals(FastFftLength.nextFiveSmooth(8319), 8640)
