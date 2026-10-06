package signal4s.fft

import gale.linalg.Vec
import signal4s.*

class SpectralBoundarySuite extends munit.FunSuite:
  test("singleton named windows have unit taps and gains in both conventions"):
    for conv <- List(WindowConvention.Periodic, WindowConvention.Symmetric);
        spec <- List(WindowSpec.Hann(1, conv), WindowSpec.Hamming(1, conv), WindowSpec.Blackman(1, conv), WindowSpec.Kaiser(1, 5.0, conv)) do
      val window = Window.fromSpec(spec).orThrow
      assertEquals(window.taps.toSeq, Seq(1.0))
      assertEquals(window.coherentGain, 1.0)
      assertEquals(window.powerGain, 1.0)
      assertEquals(window.enbwBins, 1.0)

  test("custom windows reject nonfinite taps and unrepresentable aggregate gains"):
    for tap <- List(Double.NaN, Double.PositiveInfinity, Double.NegativeInfinity, Double.MaxValue) do
      assert(Window.fromTaps(Vec(tap), WindowConvention.Periodic).isLeft)
    val large = Window.fromTaps(Vec(1e154), WindowConvention.Periodic).orThrow
    assertEqualsDouble(large.enbwBins, 1.0, 1e-12)

  test("Welch refuses undefined selected normalization but permits zero-coherent density windows"):
    val fs = SampleRate.hertz(8.0).orThrow
    val zero = Window.fromTaps(Vec.zeros(4), WindowConvention.Periodic).orThrow
    for scaling <- List(SpectralScaling.Density, SpectralScaling.Spectrum) do
      assert(WelchPlan(zero, 2, fs, 4, scaling = scaling).isLeft)
    val signed = Window.fromTaps(Vec(1.0, -1.0, 1.0, -1.0), WindowConvention.Periodic).orThrow
    assert(WelchPlan(signed, 2, fs, 4, scaling = SpectralScaling.Spectrum).isLeft)
    assert(WelchPlan(signed, 2, fs, 4, scaling = SpectralScaling.Density).isRight)
    val hugeFs = SampleRate.hertz(Double.MaxValue).orThrow
    val rectangular = Window.fromTaps(Vec(1.0, 1.0), WindowConvention.Periodic).orThrow
    assert(WelchPlan(rectangular, 2, hugeFs, 2).isLeft)

  test("singleton periodogram remains finite with either window convention"):
    val fs = SampleRate.hertz(8.0).orThrow
    for conv <- List(WindowConvention.Periodic, WindowConvention.Symmetric) do
      val window = Window.fromSpec(WindowSpec.Hann(1, conv)).orThrow
      val result = Periodogram(Vec(2.0), fs, window, 1, detrend = Detrend.None).orThrow
      assertEquals(result.power.toSeq, Seq(0.5))
