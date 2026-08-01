package signal4s.fft

import gale.linalg.{DVec, Vec}
import signal4s.*

class SpectralEdgeSuite extends munit.FunSuite:

  test("detrending preserves none, removes means, and removes linear trends"):
    val frame = Vec(1.0, 3.0, 5.0, 7.0)
    assert(DetrendOps(frame, Detrend.None) eq frame)
    assertClose(DetrendOps(frame, Detrend.Mean), Vec(-3.0, -1.0, 1.0, 3.0))
    assertClose(DetrendOps(frame, Detrend.Linear), Vec(0.0, 0.0, 0.0, 0.0))
    assertClose(DetrendOps(Vec(2.5), Detrend.Linear), Vec(0.0))

  test("Welch supports padded median aggregation and records segment metadata"):
    val fs = SampleRate.hertz(20.0).orThrow
    val window = Window.fromSpec(WindowSpec.Rectangular(4, WindowConvention.Periodic)).orThrow
    val plan = WelchPlan(
      window,
      hop = 4,
      sampleRate = fs,
      nfft = 8,
      detrend = Detrend.Linear,
      scaling = SpectralScaling.Spectrum,
      average = AverageMethod.Median
    ).orThrow
    val signal = Vec.tabulate(12)(i => 1.0 + 0.2 * i + math.sin(0.7 * i))
    val result = plan.estimate(signal).orThrow
    assertEquals(result.segmentCount, 3)
    assertEquals(result.power.length, 5)
    assertEquals(result.average, AverageMethod.Median)
    assertEquals(result.scaling, SpectralScaling.Spectrum)
    assert(result.diagnostics.contains("hop=4"))
    var i = 0
    while i < result.power.length do
      assert(result.power(i).isFinite && result.power(i) >= 0.0, s"power@$i")
      i += 1

  test("Welch rejects impossible plans, short signals, and twosided estimates"):
    val fs = SampleRate.hertz(10.0).orThrow
    val window = Window.fromSpec(WindowSpec.Hann(4, WindowConvention.Periodic)).orThrow
    assert(WelchPlan(window, hop = 2, sampleRate = fs, nfft = 2).isLeft)
    assert(WelchPlan.hann(8, fs, noverlap = Some(8)).isLeft)

    val oneSided = WelchPlan(window, hop = 2, sampleRate = fs, nfft = 4).orThrow
    assert(oneSided.estimate(Vec(1.0, 2.0, 3.0)).isLeft)
    val twosided =
      WelchPlan(window, hop = 2, sampleRate = fs, nfft = 4, sides = SpectralSides.Twosided).orThrow
    assert(twosided.estimate(Vec.tabulate(8)(_.toDouble)).isLeft)

  private def assertClose(actual: DVec, expected: DVec): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assertEqualsDouble(actual(i), expected(i), 1e-12, clue = s"i=$i")
      i += 1
