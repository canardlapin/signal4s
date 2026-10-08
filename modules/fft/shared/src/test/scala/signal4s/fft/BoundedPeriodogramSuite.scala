package signal4s.fft

import gale.linalg.Vec
import signal4s.*

class BoundedPeriodogramSuite extends munit.FunSuite:
  private val fs = SampleRate.hertz(8).orThrow
  private def window(n: Int) = Window.fromSpec(WindowSpec.Hann(n, WindowConvention.Periodic)).orThrow

  for n <- Vector(1, 2, 4, 5, 6, 17, 32) do
    test(s"owned portable length $n agrees with ordinary Welch and Parseval"):
      val x = Vec.tabulate(n)(i => math.sin(i * 0.7) + i * 0.1)
      val w = window(n)
      for d <- Vector(Detrend.None, Detrend.Mean, Detrend.Linear); scaling <- Vector(SpectralScaling.Density, SpectralScaling.Spectrum) do
        val plan = BoundedPeriodogramPlan(w, fs, n, d, scaling).orThrow
        val actual = plan.estimate(x).orThrow
        val expected = WelchPlan(w, n, fs, n, detrend = d, scaling = scaling).orThrow.estimate(x).orThrow
        for k <- 0 until actual.power.length do assertEqualsDouble(actual.power(k), expected.power(k), 1e-11)
        assertEquals(actual.degreesOfFreedom, None)
      val actual = BoundedPeriodogramPlan(w, fs, n, Detrend.None).orThrow.estimate(x).orThrow
      val energy = (0 until n).map(i => x(i) * x(i) * w.taps(i) * w.taps(i)).sum / w.taps.toSeq.map(v => v * v).sum
      assertEqualsDouble(actual.power.toSeq.sum * 8 / n, energy, 1e-11)

  test("same-shape foreign workspace and wrong input shape refuse before reuse"):
    val a = BoundedPeriodogramPlan(window(4), fs, 4).orThrow
    val b = BoundedPeriodogramPlan(window(4), fs, 4).orThrow
    val ws = a.newWorkspace()
    val retained = a.estimateInto(Vec(0, 1, 0, -1), ws).orThrow
    assert(b.estimateInto(Vec.zeros(4), ws).isLeft)
    assert(a.estimateInto(Vec.zeros(3), ws).isLeft)
    a.estimateInto(Vec.zeros(4), ws).orThrow
    assert(retained.power.toSeq.exists(_ > 0))

  test("workspace arrays and tables exactly match declared primitive capacities"):
    for n <- Vector(1, 4, 5, 17, 128) do
      val plan = BoundedPeriodogramPlan(window(n), fs, n).orThrow
      val ws = plan.newWorkspace()
      val doubles = Vector(ws.re, ws.im, ws.chirpRe, ws.chirpIm, ws.kernelRe, ws.kernelIm).map(_.length.toLong).sum +
        (ws.tables.fwdStages.toVector ++ ws.tables.invStages.toVector).map(s => s.wRe.length.toLong + s.wIm.length).sum
      assertEquals(plan.resources.workspaceBytes, BigInt(doubles) * 8 + BigInt(ws.tables.bitrev.length) * 4)

  test("nonfinite input and genuinely unrepresentable normalized power refuse"):
    val w = Window.fromSpec(WindowSpec.Rectangular(4, WindowConvention.Periodic)).orThrow
    val bounded = BoundedPeriodogramPlan(w, fs, 4, Detrend.None).orThrow
    val legacy = WelchPlan(w, 4, fs, 4, detrend = Detrend.None).orThrow
    for bad <- Vector(Double.NaN, Double.PositiveInfinity, Double.NegativeInfinity) do
      assert(bounded.estimate(Vec(0, bad, 0, 1)).isLeft)
      assert(legacy.estimate(Vec(0, bad, 0, 1)).isLeft)
    assert(bounded.estimate(Vec.fill(4)(Double.MaxValue)).isLeft)
    assert(legacy.estimate(Vec.fill(4)(Double.MaxValue)).isLeft)

  test("finite periodogram and mean survive overflowing raw magnitudes and totals"):
    val w = Window.fromSpec(WindowSpec.Rectangular(4, WindowConvention.Periodic)).orThrow
    val rate = SampleRate.hertz(4).orThrow
    val bounded = BoundedPeriodogramPlan(w, rate, 4, Detrend.None).orThrow
    assertEqualsDouble(bounded.estimate(Vec.fill(4)(1e154)).orThrow.power(0) / 1e308, 1, 1e-14)
    assertEqualsDouble(WelchPlan(w, 4, rate, 4, detrend = Detrend.None).orThrow.estimate(Vec.fill(12)(1e154)).orThrow.power(0) / 1e308, 1, 1e-14)
    val zero = BoundedPeriodogramPlan(w, rate, 4, Detrend.Mean).orThrow.estimate(Vec.fill(4)(Double.MaxValue)).orThrow
    assertEquals(zero.power.toSeq, Seq(0.0, 0.0, 0.0))

  test("odd padded last bin doubles and spectrum versus density use fixed divisors"):
    val w = Window.fromSpec(WindowSpec.Rectangular(3, WindowConvention.Periodic)).orThrow
    for n <- Vector(5, 8) do
      val d = BoundedPeriodogramPlan(w, fs, n, Detrend.None).orThrow.estimate(Vec(1, -1, 0)).orThrow
      val s = BoundedPeriodogramPlan(w, fs, n, Detrend.None, SpectralScaling.Spectrum).orThrow.estimate(Vec(1, -1, 0)).orThrow
      for k <- 0 until d.power.length do assertEqualsDouble(s.power(k), d.power(k) * 8 / 3, 1e-12)
      assertEqualsDouble(d.power.toSeq.sum * 8 / n, 2.0 / 3, 1e-12)

  test("capacity inspection admits no allocation and rejects overflowing FFT shapes"):
    val w = window(4)
    assert(BoundedPeriodogramPlan(w, fs, 3).isLeft)
    assert(BoundedPeriodogramPlan(w, fs, Int.MaxValue).isLeft)
    val large = BoundedPeriodogramPlan(w, fs, BoundedPeriodogramPlan.MaxFftLength - 1).orThrow
    assert(large.resources.workspaceBytes > BigInt(Int.MaxValue))

  test("large DC offsets preserve small centered power without rounding the absolute mean"):
    val w = Window.fromSpec(WindowSpec.Rectangular(4, WindowConvention.Periodic)).orThrow
    val plan = BoundedPeriodogramPlan(w, fs, 4, Detrend.Mean).orThrow
    val x = Vec(1e16, 1e16 + 2, 1e16 + 4, 1e16 + 6)
    val power = plan.estimate(x).orThrow.power
    assertEqualsDouble(power.toSeq.sum * 2, 5, 1e-12)
    assertEqualsDouble(power(0), 0, 1e-14)

  test("finite even raw median does not overflow its two central values"):
    val w = Window.fromSpec(WindowSpec.Rectangular(4, WindowConvention.Periodic)).orThrow
    val rate = SampleRate.hertz(4).orThrow
    val power = WelchPlan(w, 4, rate, 4, detrend = Detrend.None, average = AverageMethod.Median).orThrow
      .estimate(Vec.fill(8)(1e154)).orThrow.power
    assertEqualsDouble(power(0) / 1e308, 1, 1e-14)
