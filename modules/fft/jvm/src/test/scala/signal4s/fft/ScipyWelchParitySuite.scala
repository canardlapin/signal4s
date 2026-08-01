package signal4s.fft

import gale.linalg.DVec
import signal4s.*
import signal4s.testkit.SciPyFixture

class ScipyWelchParitySuite extends munit.FunSuite:

  test("welch density no detrend"):
    val fixture = SciPyFixture.load("smoke.welch_density_nodetrend")
    val signal = DVec.fromSeq(fixture.array("inputs", "signal"))
    val expected = DVec.fromSeq(fixture.array("expected", "power"))
    val fs = SampleRate.hertz(1000.0).orThrow
    val win = Window.fromSpec(WindowSpec.Hann(32, WindowConvention.Periodic)).orThrow
    val plan = WelchPlan(
      win,
      hop = 16,
      sampleRate = fs,
      nfft = 32,
      detrend = Detrend.None,
      scaling = SpectralScaling.Density
    ).orThrow
    val result = plan.estimate(signal).orThrow
    assertEquals(result.degreesOfFreedom, None)
    assertClose(result.power, expected, fixture.rtol, fixture.atol)

  test("welch spectrum detrend mean"):
    val fixture = SciPyFixture.load("smoke.welch_spectrum_detrend_mean")
    val signal = DVec.fromSeq(fixture.array("inputs", "signal"))
    val expected = DVec.fromSeq(fixture.array("expected", "power"))
    val fs = SampleRate.hertz(1000.0).orThrow
    val win = Window.fromSpec(WindowSpec.Hann(32, WindowConvention.Periodic)).orThrow
    val plan = WelchPlan(
      win,
      hop = 16,
      sampleRate = fs,
      nfft = 32,
      detrend = Detrend.Mean,
      scaling = SpectralScaling.Spectrum
    ).orThrow
    val result = plan.estimate(signal).orThrow
    assertClose(result.power, expected, fixture.rtol, fixture.atol)

  test("periodogram spectrum hann"):
    val fixture = SciPyFixture.load("smoke.periodogram_spectrum_hann")
    val signal = DVec.fromSeq(fixture.array("inputs", "signal"))
    val expected = DVec.fromSeq(fixture.array("expected", "power"))
    val fs = SampleRate.hertz(1000.0).orThrow
    val result =
      Periodogram
        .hann(signal, fs, nfft = Some(128), detrend = Detrend.None, scaling = SpectralScaling.Spectrum)
        .orThrow
    assertClose(result.power, expected, fixture.rtol, fixture.atol)

  private def assertClose(a: DVec, b: DVec, rtol: Double, atol: Double): Unit =
    assertEquals(a.length, b.length)
    var i = 0
    while i < a.length do
      val tol = atol + rtol * math.max(math.abs(a(i)), math.abs(b(i)))
      assert(math.abs(a(i) - b(i)) <= tol, s"i=$i ${a(i)} vs ${b(i)}")
      i += 1
