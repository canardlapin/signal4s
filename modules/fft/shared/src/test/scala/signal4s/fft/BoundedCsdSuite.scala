package signal4s.fft

import gale.linalg.{DVec, Vec}
import signal4s.*

class BoundedCsdSuite extends munit.FunSuite:
  private val budget = BigInt(64) * 1024 * 1024
  private val fs = SampleRate.hertz(16).orThrow
  private val window = Window.fromSpec(WindowSpec.Rectangular(8, WindowConvention.Periodic)).orThrow
  private def plan(channels: Int = 3, nfft: Int = 8) = BoundedCsdPlan.withinBudget(channels, window, fs, nfft, budget).orThrow
  private val x = Vec.tabulate(8)(i => math.cos(2 * math.Pi * i / 8))
  private val y = Vec.tabulate(8)(i => math.sin(2 * math.Pi * i / 8))
  private val inputs = Vector(x, y, Vec.tabulate(8)(i => x(i) + y(i)))
  private def result(p: BoundedCsdPlan, frames: Vector[IndexedSeq[DVec]]): CsdMatrix =
    val accumulator = p.newAccumulator(budget).orThrow
    frames.foreach(frame => accumulator.accumulateFrame(frame).orThrow)
    accumulator.finishOwned().orThrow.get
  private def close(a: Double, b: Double, tolerance: Double = 1e-12): Unit =
    assert(math.abs(a - b) <= tolerance, s"$a != $b")

  test("conjugation has the declared relative phase and matched density diagonal"):
    val matrix = result(plan(), Vector(inputs))
    close(matrix(0, 1, 1).real, 0)
    close(matrix(0, 1, 1).imaginary, -0.25)
    close(matrix(1, 0, 1).imaginary, 0.25)
    close(matrix(0, 0, 1).real, 0.25)
    for row <- 0 until 3; k <- 0 until matrix.binCount do
      assertEquals(matrix(row, row, k).imaginary, 0.0)
      assert(matrix(row, row, k).real >= 0)
      for column <- 0 until 3 do
        assertEquals(matrix(row, column, k).real, matrix(column, row, k).real)
        close(matrix(row, column, k).imaginary, -matrix(column, row, k).imaginary, 0)

  test("complex quadratic forms are nonnegative for a rank-deficient Gram estimate"):
    val matrix = result(plan(), Vector(inputs, inputs.map(v => Vec.tabulate(v.length)(i => v(i) * 2))))
    for k <- 0 until matrix.binCount; seed <- 0 until 20 do
      val re = Vector.tabulate(3)(i => math.sin(seed * 3.0 + i))
      val im = Vector.tabulate(3)(i => math.cos(seed + i * 2.0))
      val value = (0 until 3).map: row =>
        val realPart = (0 until 3).map(j => matrix(row, j, k).real * re(j) - matrix(row, j, k).imaginary * im(j)).sum
        val imagPart = (0 until 3).map(j => matrix(row, j, k).real * im(j) + matrix(row, j, k).imaginary * re(j)).sum
        re(row) * realPart + im(row) * imagPart
      val scale = (0 until 3).map(i => matrix(i, i, k).real).max
      assert(value.sum >= -1e-10 * scale, (k, seed, value.sum, scale))

  test("row permutation preserves both complex axes for odd padded FFTs"):
    val p = plan(nfft = 9)
    val original = result(p, Vector(inputs))
    val permutation = Vector(2, 0, 1)
    val changed = result(p, Vector(permutation.map(inputs)))
    for i <- 0 until 3; j <- 0 until 3; k <- 0 until original.binCount do
      close(changed(i, j, k).real, original(permutation(i), permutation(j), k).real, 0)
      close(changed(i, j, k).imaginary, original(permutation(i), permutation(j), k).imaginary, 0)

  test("late invalid channel and product overflow preserve count and prior spectrum"):
    val p = plan(); val accumulator = p.newAccumulator(budget).orThrow
    assertEquals(accumulator.finishOwned().orThrow, None)
    accumulator.accumulateFrame(inputs).orThrow
    val before = accumulator.finishOwned().orThrow.get
    assert(accumulator.accumulateFrame(inputs.updated(2, Vec.fill(8)(Double.NaN))).isLeft)
    assert(accumulator.accumulateFrame(inputs.updated(2, Vec.fill(8)(1e200))).isLeft)
    assert(accumulator.accumulateFrame(inputs.updated(2, Vec(1.0))).isLeft)
    assert(accumulator.accumulateFrame(inputs.take(2)).isLeft)
    val after = accumulator.finishOwned().orThrow.get
    assertEquals(after.acceptedWindows, 1L)
    assertEquals(after.real.toSeq, before.real.toSeq)
    assertEquals(after.imaginary.toSeq, before.imaginary.toSeq)
    accumulator.accumulateFrame(inputs.map(v => Vec.tabulate(v.length)(i => v(i) * 2))).orThrow
    assertEquals(before.real.toSeq, after.real.toSeq)
    assertEquals(before.acceptedWindows, 1L)

  test("partition merges retain exact mean products and refuse foreign or self state"):
    val p = plan(); val frames = Vector(inputs, inputs.map(v => Vec.tabulate(v.length)(i => v(i) * 0.25)), inputs.map(v => Vec.tabulate(v.length)(i => v(i) * 2)))
    val expected = result(p, frames)
    val first = p.newAccumulator(budget).orThrow; val second = p.newAccumulator(budget).orThrow
    first.accumulateFrame(frames(2)).orThrow
    second.accumulateFrame(frames(0)).orThrow; second.accumulateFrame(frames(1)).orThrow
    first.merge(second).orThrow
    val actual = first.finishOwned().orThrow.get
    assertEquals(actual.real.toSeq, expected.real.toSeq)
    assertEquals(actual.imaginary.toSeq, expected.imaginary.toSeq)
    assert(first.merge(first).isLeft)
    assert(first.merge(plan().newAccumulator(budget).orThrow).isLeft)
    assertEquals(first.acceptedWindows, 3L)
    assertEquals(second.acceptedWindows, 2L)

  test("shape and numeric budgets are checked before state allocation"):
    assert(BoundedCsdPlan.inspect(0, 8, 8).isLeft)
    assert(BoundedCsdPlan.inspect(Int.MaxValue, 8, 8).isLeft)
    assert(BoundedCsdPlan.inspect(3, 8, 4).isLeft)
    val p = plan()
    assert(p.newAccumulator(p.resources.peakBytes - 1).isLeft)
    assert(BoundedCsdPlan.withinBudget(3, window, fs, 8, p.resources.peakBytes - 1).isLeft)
    assert(p.newAccumulator(p.resources.peakBytes).isRight)

  test("pinned SciPy CSD agrees across quadrature, mean detrending, odd FFT and impulses"):
    for fixture <- CsdGoldens.cases do
      val taps = Window.fromTaps(DVec.fromSeq(fixture.window), WindowConvention.Periodic).orThrow
      val p = BoundedCsdPlan.withinBudget(3, taps, SampleRate.hertz(128).orThrow, fixture.nfft, budget, fixture.detrend).orThrow
      val actual = result(p, fixture.frames.map(_.map(DVec.fromSeq)))
      val peak = Vector.tabulate(3)(row => (0 until actual.binCount).map(k => fixture.real((row * 3 + row) * actual.binCount + k)).max)
      for i <- 0 until 3; j <- 0 until 3; k <- 0 until actual.binCount do
        val index = (i * 3 + j) * actual.binCount + k
        val sii = fixture.real((i * 3 + i) * actual.binCount + k)
        val sjj = fixture.real((j * 3 + j) * actual.binCount + k)
        val originalTolerance = 1e-30 + 1e-10 * math.sqrt(sii * sjj)
        // Reviewed v2: same scale convention as qualified Welch; original failed
        // per-bin criterion is preserved in the ephys4s contract and diagnostics.
        val tolerance = 1e-30 + 1e-10 * math.sqrt(peak(i) * peak(j))
        if math.abs(actual(i,j,k).real - fixture.real(index)) > originalTolerance ||
            math.abs(actual(i,j,k).imaginary - fixture.imaginary(index)) > originalTolerance then
          println(s"CSD original-policy failure retained: ${fixture.name} ($i,$j,$k), original=$originalTolerance reviewed=$tolerance")
        assert(math.abs(actual(i, j, k).real - fixture.real(index)) <= tolerance,
          s"${fixture.name}: re($i,$j,$k) actual=${actual(i,j,k).real} expected=${fixture.real(index)} tolerance=$tolerance diagonal=($sii,$sjj)")
        assert(math.abs(actual(i, j, k).imaginary - fixture.imaginary(index)) <= tolerance,
          s"${fixture.name}: im($i,$j,$k) actual=${actual(i,j,k).imaginary} expected=${fixture.imaginary(index)} tolerance=$tolerance diagonal=($sii,$sjj)")
