package signal4s.fft

import gale.linalg.{DVec, Vec}

class LogPowerBaselineSuite extends munit.FunSuite:
  private def ok[A](value: Either[?, A]): A = value.fold(e => fail(e.toString), identity)
  private def baseline(convention: LogPowerConvention, values: Vector[Double] = Vector(1.0, 4.0),
      weights: Vector[Double] = Vector(1.0, 1.0), scale: LogPowerScale = LogPowerScale.Decibels): LogPowerBaseline =
    ok(LogPowerBaseline.fit(DVec.fromSeq(values), DVec.fromSeq(weights), Vector.fill(values.size)(true), convention, scale))

  test("individual matched targets retain the Jensen gap and mean-log zero"):
    val legacy = baseline(LogPowerConvention.LegacyLogRatio); val meanLog = baseline(LogPowerConvention.MeanLogRatio)
    val targets = Vector(1.0, 4.0)
    assertEqualsDouble(targets.map(p => ok(legacy.normalize(p))).sum / 2, -0.969100130080564, 2e-14)
    assertEqualsDouble(targets.map(p => ok(meanLog.normalize(p))).sum / 2, 0.0, 2e-14)
    assertEqualsDouble(ok(legacy.normalize(2.5)), 0.0, 0.0)
    assertEqualsDouble(ok(meanLog.normalize(2.5)), 0.969100130080564, 2e-14)

  test("ten percent power increase preserves both explicitly different estimands"):
    val targets = Vector(1.1, 4.4)
    assertEqualsDouble(targets.map(p => ok(baseline(LogPowerConvention.LegacyLogRatio).normalize(p))).sum / 2, -0.5551732784983132, 2e-14)
    assertEqualsDouble(targets.map(p => ok(baseline(LogPowerConvention.MeanLogRatio).normalize(p))).sum / 2, 0.4139268515822505, 2e-14)

  test("nonuniform weights agree with independently computed arithmetic and geometric baselines"):
    val legacy = baseline(LogPowerConvention.LegacyLogRatio, weights = Vector(1.0, 3.0))
    val meanLog = baseline(LogPowerConvention.MeanLogRatio, weights = Vector(1.0, 3.0))
    assertEqualsDouble(ok(legacy.normalize(2)), 10 * math.log10(2.0 / 3.25), 2e-14)
    assertEqualsDouble(ok(meanLog.normalize(2)), -5 * math.log10(2), 2e-14)
    assertEquals(meanLog.selectedCount, 2); assertEquals(meanLog.positiveWeightCount, 2); assertEquals(meanLog.weightTotal, 4.0)

  test("scale and base declare dimensionless natural log log10 and power dB"):
    for scale <- LogPowerScale.values do
      val b = baseline(LogPowerConvention.MeanLogRatio, Vector(1.0, 1.0), scale = scale)
      val expected = if scale == LogPowerScale.NaturalLog then math.log(10) else scale.multiplier
      assertEqualsDouble(ok(b.normalize(10)), expected, 2e-14)
      assertEquals(scale.base, if scale == LogPowerScale.NaturalLog then math.E else 10.0)

  test("masked absence and zero weights retain distinct observation counts"):
    val b = ok(LogPowerBaseline.fit(Vec(1, 4, 9, Double.NaN), Vec(1, 3, 0, 1), Vector(true, true, true, false),
      LogPowerConvention.MeanLogRatio, LogPowerScale.Decibels))
    assertEquals(b.selectedCount, 3); assertEquals(b.positiveWeightCount, 2); assertEquals(b.weightTotal, 4.0)
    assertEqualsDouble(ok(b.normalize(2)), -5 * math.log10(2), 2e-14)

  test("nonpositive nonfinite empty and all-omitted baselines refuse without flooring"):
    for value <- Vector(0.0, -1.0, Double.NaN, Double.PositiveInfinity) do
      assert(LogPowerBaseline.fit(Vec(value), Vec(1), Vector(true), LogPowerConvention.MeanLogRatio, LogPowerScale.Decibels).isLeft)
      assert(baseline(LogPowerConvention.MeanLogRatio).normalize(value).isLeft)
    for mask <- Vector(Vector(true), Vector(false)) do
      assert(LogPowerBaseline.fit(Vec(1), Vec(0), mask, LogPowerConvention.MeanLogRatio, LogPowerScale.Decibels).isLeft)
    assert(LogPowerBaseline.fit(Vec.zeros(0), Vec.zeros(0), Vector.empty, LogPowerConvention.MeanLogRatio, LogPowerScale.Decibels).isLeft)

  test("weight length validity and arithmetic capacity are explicit admission checks"):
    for weight <- Vector(-1.0, Double.NaN, Double.PositiveInfinity) do
      assert(LogPowerBaseline.fit(Vec(1), Vec(weight), Vector(true), LogPowerConvention.MeanLogRatio, LogPowerScale.Decibels).isLeft)
    assert(LogPowerBaseline.fit(Vec(1), Vec(1, 2), Vector(true), LogPowerConvention.MeanLogRatio, LogPowerScale.Decibels).isLeft)
    assert(LogPowerBaseline.fit(Vec(Double.MaxValue), Vec(2), Vector(true), LogPowerConvention.LegacyLogRatio, LogPowerScale.Decibels).isLeft)
    val large = baseline(LogPowerConvention.LegacyLogRatio, Vector(Double.MaxValue, Double.MaxValue))
    assertEqualsDouble(ok(large.normalize(Double.MaxValue)), 0, 0)

  test("permutations and repeated equal weighted entries preserve the fitted baseline"):
    for convention <- LogPowerConvention.values do
      val original = baseline(convention, Vector(1, 4, 3), Vector(1, 2, 3))
      val reversed = baseline(convention, Vector(3, 4, 1), Vector(3, 2, 1))
      assertEquals(original.baselineLogPower, reversed.baselineLogPower)
      val held = original.baselineLogPower
      baseline(convention, Vector(1, 9), Vector(1, 1))
      assertEquals(original.baselineLogPower, held)

  test("matched log distributions with different variability retain the declared contrast"):
    val tight = baseline(LogPowerConvention.MeanLogRatio, Vector(1, 1))
    val variable = baseline(LogPowerConvention.MeanLogRatio, Vector(0.25, 4))
    assertEqualsDouble(tight.baselineLogPower, variable.baselineLogPower, 2e-14)
    val legacy = baseline(LogPowerConvention.LegacyLogRatio, Vector(0.25, 4))
    assert(ok(legacy.normalize(1)) < -3)
    assertEqualsDouble(ok(variable.normalize(1)), 0, 2e-14)
