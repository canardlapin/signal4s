package signal4s.fft

import gale.linalg.{DVec, Vec}
import signal4s.*

/** Cost-model choices must stay numerically equivalent to direct convolution. */
class AutoConvolutionPlanSuite extends munit.FunSuite:

  FftBackend.ensureInstalled()

  test("Auto selects FFT for a realistic dense FIR and matches direct"):
    val signal = Vec.tabulate(1024)(i => math.sin(0.031 * i) + 0.02 * (i % 11))
    val kernel = Kernel.causal(Vec.tabulate(512)(i => math.cos(0.017 * i) / (i + 1.0))).orThrow
    val plan = Convolution.plan(kernel, signal.length, OutputRegion.Full, ConvolutionMethod.Auto).orThrow
    assertEquals(plan.selectedMethod, ConvolutionMethod.Fft)
    assertClose(
      plan(signal).orThrow,
      Convolution(signal, kernel, OutputRegion.Full, ConvolutionMethod.Direct).orThrow
    )

  test("automatic overlap-add block handles a ragged final block"):
    val signal = Vec.tabulate(137)(i => math.sin(0.19 * i) - 0.03 * (i % 4))
    val kernel = Kernel.causal(Vec.tabulate(9)(i => 1.0 / (i + 1.0))).orThrow
    val plan =
      Convolution.plan(kernel, signal.length, OutputRegion.Full, ConvolutionMethod.OverlapAdd(0)).orThrow
    assert(
      plan.selectedMethod match
        case ConvolutionMethod.OverlapAdd(block) => block > 0
        case _                                    => false,
      clue = s"selected ${plan.selectedMethod}"
    )
    assertClose(
      plan(signal).orThrow,
      Convolution(signal, kernel, OutputRegion.Full, ConvolutionMethod.Direct).orThrow
    )

  private def assertClose(actual: DVec, expected: DVec): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assertEqualsDouble(actual(i), expected(i), 1e-8, clue = s"i=$i")
      i += 1
