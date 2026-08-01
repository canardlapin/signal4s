package signal4s.fft

import gale.linalg.{DVec, Vec}
import signal4s.*

/** STFT construction and synthesis validation must fail before numeric work begins. */
class StftContractSuite extends munit.FunSuite:

  private val fs = SampleRate.hertz(16.0).orThrow
  private val window = Window.fromSpec(WindowSpec.Rectangular(4, WindowConvention.Periodic)).orThrow

  test("canonical dual has the expected overlap normalization and validates hop"):
    assert(CanonicalDual.fromAnalysis(window, hop = 0).isLeft)
    val dual = CanonicalDual.fromAnalysis(window, hop = 2).orThrow
    assertEquals(dual.length, 4)
    assertEquals(dual.hop, 2)
    var i = 0
    while i < dual.length do
      assertEqualsDouble(dual.taps(i), 0.5, 1e-15, clue = s"tap=$i")
      i += 1

  test("STFT construction validates FFT and custom synthesis dimensions"):
    assert(StftPlan(window, hop = 2, sampleRate = fs, nfft = 2).isLeft)
    val shortSynthesis = Window.fromSpec(WindowSpec.Rectangular(3, WindowConvention.Periodic)).orThrow
    assert(StftPlan(window, hop = 2, sampleRate = fs, nfft = 4, synthesis = Some(shortSynthesis)).isLeft)
    val plan = StftPlan(window, hop = 2, sampleRate = fs, nfft = 4, synthesis = Some(window)).orThrow
    assertEquals(plan.dual, None)
    assertEquals(plan.hop, 2)
    assertEquals(plan.frameLength, 4)

  test("analysis and synthesis reject empty, incompatible, twosided, and nonpositive shapes"):
    val plan = StftPlan(window, hop = 2, sampleRate = fs, nfft = 4).orThrow
    assert(plan.analyze(DVec.zeros(0)).isLeft)
    val tf = plan.analyze(Vec.tabulate(8)(i => math.sin(0.4 * i))).orThrow
    assert(tf.frameCount > 0)
    assertEquals(tf.binCount, 3)

    val wrongNfft = TimeFrequency(tf.columns, tf.times, tf.frequencies, nfft = 8, onesided = true).orThrow
    assert(plan.synthesize(wrongNfft, outputLength = 8).isLeft)
    val twosided = TimeFrequency(tf.columns, tf.times, tf.frequencies, nfft = 4, onesided = false).orThrow
    assert(plan.synthesize(twosided, outputLength = 8).isLeft)
    assert(plan.synthesize(tf, outputLength = 0).isLeft)
