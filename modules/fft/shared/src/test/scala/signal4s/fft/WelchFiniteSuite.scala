package signal4s.fft

import gale.linalg.Vec
import signal4s.*

class WelchFiniteSuite extends munit.FunSuite:
  private val window = Window.fromSpec(WindowSpec.Rectangular(4, WindowConvention.Periodic)).orThrow
  private val fs = SampleRate.hertz(4).orThrow

  test("finite normalized density survives overflowing raw squared FFT magnitude"):
    val result = WelchPlan(window, 4, fs, 4, detrend = Detrend.None).orThrow.estimate(Vec.fill(4)(1e154)).orThrow
    assert(result.power(0).isFinite)
    assertEqualsDouble(result.power(0) / 1e308, 1.0, 1e-14)
    assertEquals(result.power.toSeq.drop(1), Seq(0.0, 0.0))

  test("constant maximum finite values detrend to zero without raw sum overflow"):
    val result = WelchPlan(window, 4, fs, 4, detrend = Detrend.Mean).orThrow.estimate(Vec.fill(4)(Double.MaxValue)).orThrow
    assertEquals(result.power.toSeq, Seq(0.0, 0.0, 0.0))
