package signal4s.fft

import gale.linalg.{DVec, Vec}
import signal4s.*

class BoundedStftSuite extends munit.FunSuite:
  private val fs = SampleRate.hertz(16.0).orThrow
  private def rect(n: Int) = Window.fromSpec(WindowSpec.Rectangular(n, WindowConvention.Periodic)).orThrow
  private def transform(plan: BoundedStftPlan, input: DVec) = plan.analyze(input, plan.resources.peakBytes).orThrow
  private def close(actual: Double, expected: Double, tolerance: Double = 1e-12): Unit =
    assert(math.abs(actual - expected) <= tolerance * math.max(1.0, math.abs(expected)), s"$actual != $expected")

  test("analysis-only Hann accepts the noninvertible no-overlap window"):
    val window = Window.fromSpec(WindowSpec.Hann(4, WindowConvention.Periodic)).orThrow
    assert(StftPlan(window, 4, fs, 4).isLeft)
    val plan = BoundedStftPlan(window, fs, 4).orThrow
    val frame = transform(plan, Vec(0.0, 0.0, 1.0, 0.0))
    assertEquals(frame.bins.real.toSeq, Seq(1.0, -1.0, 1.0))
    close(frame.bins.imaginary(1), 0.0)

  test("frame-start impulse phase is preserved without a centering rotation"):
    val frame = transform(BoundedStftPlan(rect(4), fs, 4).orThrow, Vec(0.0, 1.0, 0.0, 0.0))
    close(frame.bins.real(0), 1.0); close(frame.bins.imaginary(0), 0.0)
    close(frame.bins.real(1), 0.0); close(frame.bins.imaginary(1), -1.0)
    close(frame.bins.real(2), -1.0); close(frame.bins.imaginary(2), 0.0)
    assertEquals(frame.phaseOrigin, "frame_start; raw rFFT; phase_shift=None")

  test("complex magnitude and density do not silently double interior bins"):
    val input = Vec.tabulate(8)(i => math.cos(2.0 * math.Pi * i / 8))
    val magnitude = transform(BoundedStftPlan(rect(8), fs, 8, scaling = StftScaling.Magnitude).orThrow, input)
    val density = transform(BoundedStftPlan(rect(8), fs, 8, scaling = StftScaling.Density).orThrow, input)
    close(magnitude.bins.real(1), 0.5)
    close(density.bins.real(1), 4.0 / math.sqrt(128.0))

  test("raw folded Parseval holds for odd prime and padded FFT lengths"):
    for nfft <- Vector(5, 7, 8, 11, 16) do
      val input = Vec(1.0, -2.0, 0.5, 3.0, -0.25)
      val frame = transform(BoundedStftPlan(rect(5), fs, nfft).orThrow, input)
      var spectral = 0.0
      for k <- 0 until frame.bins.length do
        val factor = if k == 0 || (nfft % 2 == 0 && k == nfft / 2) then 1.0 else 2.0
        spectral += factor * (frame.bins.real(k) * frame.bins.real(k) + frame.bins.imaginary(k) * frame.bins.imaginary(k))
      close(spectral / nfft, input.toSeq.map(x => x * x).sum)

  test("window normalization preserves signed phase and cancels window scale"):
    for scale <- Vector(1e-150, 1.0, 1e150); sign <- Vector(-1.0, 1.0) do
      val window = Window.fromTaps(Vec.tabulate(4)(_ => sign * scale), WindowConvention.Periodic).orThrow
      val input = Vec(1.0, 1.0, 1.0, 1.0)
      close(transform(BoundedStftPlan(window, fs, 4, scaling = StftScaling.Magnitude).orThrow, input).bins.real(0), sign)
      close(transform(BoundedStftPlan(window, fs, 4, scaling = StftScaling.Density).orThrow, input).bins.real(0), sign * 0.5)

  test("density restoration avoids overflowing raw amplitude and power products"):
    val plan = BoundedStftPlan(rect(4), fs, 4, scaling = StftScaling.Density).orThrow
    val actual = transform(plan, Vec.tabulate(4)(_ => Double.MaxValue / 2)).bins.real(0)
    assert(actual.isFinite)
    close(actual / (Double.MaxValue / 2), 0.5)
    assert(BoundedStftPlan(rect(4), fs, 4).orThrow.analyze(Vec.tabulate(4)(_ => Double.MaxValue), plan.resources.peakBytes).isLeft)

  test("density survives underflowed raw window energy and extreme finite sample rates"):
    val tiny = Window.fromTaps(Vec.fill(4)(1e-300), WindowConvention.Periodic).orThrow
    assertEquals(tiny.powerGain, 0.0)
    close(transform(BoundedStftPlan(tiny, fs, 4, scaling = StftScaling.Density).orThrow,
      Vec.fill(4)(1.0)).bins.real(0), 0.5)
    val slow = SampleRate.hertz(1e-300).orThrow
    val output = transform(BoundedStftPlan(rect(4), slow, 4, scaling = StftScaling.Density).orThrow,
      Vec.fill(4)(1e-150)).bins.real(0)
    close(output, 2.0)
    val fastest = SampleRate.hertz(Double.MaxValue).orThrow
    val high = transform(BoundedStftPlan(rect(4), fastest, 4, scaling = StftScaling.Density).orThrow,
      Vec.fill(4)(1e150)).bins.real(0)
    close(high * math.sqrt(Double.MaxValue) / 1e150, 2.0)
    assert(BoundedStftPlan(rect(4), SampleRate.hertz(Double.MinPositiveValue).orThrow, 4).isLeft)

  test("returned columns stay detached after workspace reuse and failures"):
    val plan = BoundedStftPlan(rect(4), fs, 7).orThrow
    val workspace = plan.newWorkspace(plan.resources.peakBytes).orThrow
    val first = plan.analyzeInto(Vec(1.0, 2.0, 3.0, 4.0), workspace).orThrow
    val snapshot = first.bins.copy
    val second = plan.analyzeInto(Vec(4.0, 3.0, 2.0, 1.0), workspace).orThrow
    assert(first.bins.real.toSeq != second.bins.real.toSeq)
    assertEquals(first.bins.real.toSeq, snapshot.real.toSeq)
    assertEquals(first.bins.imaginary.toSeq, snapshot.imaginary.toSeq)
    assert(plan.analyzeInto(Vec(1.0), workspace).isLeft)
    assert(plan.analyzeInto(Vec(1.0, Double.NaN, 2.0, 3.0), workspace).isLeft)
    val foreign = BoundedStftPlan(rect(4), fs, 7).orThrow
    assert(foreign.analyzeInto(Vec(1.0, 2.0, 3.0, 4.0), workspace).isLeft)
    assertEquals(first.bins.real.toSeq, snapshot.real.toSeq)

  test("budget and maximum shape inspection refuse before allocating FFT arrays"):
    val huge = BoundedStftPlan.inspect(4, BoundedStftPlan.MaxFftLength).orThrow
    assert(huge.workspaceBytes > BigInt(Int.MaxValue))
    assert(BoundedStftPlan.inspect(4, Int.MaxValue).isLeft)
    assert(BoundedStftPlan.inspect(0, 8).isLeft)
    assert(BoundedStftPlan.inspect(8, 4).isLeft)
    val plan = BoundedStftPlan(rect(4), fs, 7).orThrow
    assert(plan.newWorkspace(plan.resources.peakBytes - 1).isLeft)
    assert(plan.newWorkspace(plan.resources.peakBytes).isRight)
    assertEquals(plan.resources.fftWorkingLength, 16)
    assert(BoundedStftPlan.withinBudget(rect(4), fs, 7, 0).isLeft)
    val budget = plan.resources.peakBytes.max(plan.resources.constructionAdditionalBytes + plan.resources.windowBytes)
    assert(BoundedStftPlan.withinBudget(rect(4), fs, 7, budget).isRight)

  test("undefined normalization refuses; a raw zero window remains valid analysis"):
    val zero = Window.fromTaps(Vec.zeros(4), WindowConvention.Periodic).orThrow
    assert(BoundedStftPlan(zero, fs, 4, scaling = StftScaling.Magnitude).isLeft)
    assert(BoundedStftPlan(zero, fs, 4, scaling = StftScaling.Density).isLeft)
    val frame = transform(BoundedStftPlan(zero, fs, 4).orThrow, Vec(1.0, 2.0, 3.0, 4.0))
    assert(frame.bins.real.toSeq.forall(_ == 0.0))
    val signed = Window.fromTaps(Vec(1.0, -1.0, 1.0, -1.0), WindowConvention.Periodic).orThrow
    assert(BoundedStftPlan(signed, fs, 4, scaling = StftScaling.Magnitude).isLeft)
    assert(BoundedStftPlan(signed, fs, 4, scaling = StftScaling.Density).isRight)

  test("per-frame mean and linear detrending precede windowing"):
    for mode <- Vector(Detrend.Mean, Detrend.Linear) do
      val values = if mode == Detrend.Mean then Vec.fill(5)(3.0) else Vec.tabulate(5)(i => 2.0 + 3.0 * i)
      val result = transform(BoundedStftPlan(rect(5), fs, 7, detrend = mode).orThrow, values)
      assert(result.bins.real.toSeq.forall(v => math.abs(v) < 1e-12))
      assert(result.bins.imaginary.toSeq.forall(v => math.abs(v) < 1e-12))
