package signal4s.fft

import gale.linalg.Vec
import signal4s.*

/** FFT convolution must match direct across varied lengths (not just fixture sizes). */
class FftDirectParitySuite extends munit.FunSuite:

  override def munitFixtures = Nil

  private val shapes: List[(Int, Int)] =
    List(
      (1, 1),
      (7, 3),
      (16, 5),
      (31, 8),
      (64, 15),
      (100, 21),
      (127, 32),
      (256, 16),
      (511, 48),
      (1024, 33),
      (2048, 64),
      (3000, 50),
      (4096, 17),
      (5000, 100)
    )

  test("FFT Full plan matches Direct for varied (n,m)"):
    FftBackend.ensureInstalled()
    shapes.foreach { case (n, m) =>
      val x = Vec.tabulate(n)(i => math.sin(0.07 * i) + 0.01 * (i % 5))
      val k = Kernel.causal(Vec.tabulate(m)(i => 1.0 / (i + 1.5))).orThrow
      val direct = Convolution(x, k, OutputRegion.Full, ConvolutionMethod.Direct).orThrow
      val plan = Convolution.plan(k, n, OutputRegion.Full, ConvolutionMethod.Fft).orThrow
      val fft = plan(x).orThrow
      assertEquals(fft.length, direct.length, clue = s"n=$n m=$m length")
      var i = 0
      var maxAbs = 0.0
      while i < fft.length do
        maxAbs = math.max(maxAbs, math.abs(fft(i) - direct(i)))
        i += 1
      assert(
        maxAbs <= 1e-7,
        clue =
          s"n=$n m=$m nfft=${signal4s.fft.internal.FastFftLength(n + m - 1)} maxAbs=$maxAbs"
      )
    }

  test("FFT Valid and Input(Zero) match Direct"):
    FftBackend.ensureInstalled()
    List((64, 9), (200, 25), (1024, 41)).foreach { case (n, m) =>
      val x = Vec.tabulate(n)(i => math.cos(0.11 * i))
      val k = Kernel.centeredOdd(Vec.tabulate(m)(i => if i == m / 2 then 1.0 else 0.1)).orThrow
      List(
        OutputRegion.Valid,
        OutputRegion.Input(Boundary.Zero)
      ).foreach { region =>
        val direct = Convolution(x, k, region, ConvolutionMethod.Direct).orThrow
        val plan = Convolution.plan(k, n, region, ConvolutionMethod.Fft).orThrow
        val fft = plan(x).orThrow
        assertEquals(fft.length, direct.length, clue = s"$region n=$n")
        var i = 0
        while i < fft.length do
          assertEqualsDouble(fft(i), direct(i), 1e-8, clue = s"$region n=$n i=$i")
          i += 1
      }
    }
