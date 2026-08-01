package signal4s.fft

import gale.linalg.{DVec, Vec}
import signal4s.*
import signal4s.internal.DirectConvolution

class CircularFftSuite extends munit.FunSuite:

  FftBackend.ensureInstalled()

  test("FFT circular convolution retains direct semantics for causal and centered kernels"):
    val signal = Vec(1.0, -2.0, 0.5, 3.0, -1.0, 2.0, 0.25)
    val kernels = List(
      Kernel.causal(Vec(1.0, 0.5, -0.25)).orThrow,
      Kernel.centeredOdd(Vec(0.25, 0.5, 0.25)).orThrow,
      Kernel.at(Vec(1.0, -1.0, 0.5, 0.25), zeroLagIndex = 2).orThrow
    )
    kernels.foreach { kernel =>
      val expected = DirectConvolution.circular(signal, kernel, signal.length)
      val actual = FftBackend.circular(signal, kernel, signal.length).orThrow
      assertClose(actual, expected, clue = s"origin=${kernel.zeroLagIndex}")
    }

  test("packed real circular convolution matches direct over practical even periods"):
    List(8, 64, 256).foreach { period =>
      val signal = Vec.tabulate(period)(i => math.sin(0.13 * i) - 0.02 * (i % 7))
      val kernel = Kernel.at(
        Vec.tabulate(17)(i => math.cos(0.19 * i) / (i + 1.0)),
        zeroLagIndex = 8
      ).orThrow
      val expected = DirectConvolution.circular(signal, kernel, period)
      val actual = FftBackend.circular(signal, kernel, period).orThrow
      assertClose(actual, expected, clue = s"period=$period")
      assertClose(
        Convolution.circular(signal, kernel, period).orThrow,
        expected,
        clue = s"automatic period=$period"
      )
    }

  test("FFT circular convolution validates the same period contract as direct"):
    val signal = Vec(1.0, 2.0, 3.0)
    val kernel = Kernel.causal(Vec(1.0)).orThrow
    assertEquals(FftBackend.circular(signal, kernel, 0), Left(SignalError.InvalidPeriod(0)))
    assertEquals(
      FftBackend.circular(signal, kernel, 2),
      Left(SignalError.LengthMismatch(expected = 2, actual = 3))
    )

  private def assertClose(actual: DVec, expected: DVec, clue: String): Unit =
    assertEquals(actual.length, expected.length, clue)
    var i = 0
    while i < actual.length do
      assertEqualsDouble(actual(i), expected(i), 1e-10, clue = s"$clue@$i")
      i += 1
