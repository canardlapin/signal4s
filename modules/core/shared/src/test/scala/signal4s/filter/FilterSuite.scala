package signal4s.filter

import gale.linalg.{DVec, Vec}
import signal4s.*

class FilterSuite extends munit.FunSuite:

  test("Biquad rejects non-finite coefficients"):
    assert(Biquad(1.0, 0.0, 0.0, Double.NaN, 0.0).isLeft)

  test("SciPy SOS rows normalize a0 and reject malformed sections"):
    val normalized = Biquad.fromSosRow(IArray(2.0, 4.0, 6.0, 2.0, -1.0, 0.5)).orThrow
    assertEqualsDouble(normalized.b0, 1.0, 1e-15)
    assertEqualsDouble(normalized.b1, 2.0, 1e-15)
    assertEqualsDouble(normalized.a1, -0.5, 1e-15)
    assert(Biquad.fromSosRow(IArray(1.0, 0.0)).isLeft)
    assert(Biquad.fromSosRow(IArray(1.0, 0.0, 0.0, 0.0, 0.0, 0.0)).isLeft)
    assert(
      SecondOrderCascade
        .fromSosMatrix(IArray(IArray(1.0, 0.0, 0.0, 1.0, 0.0, 0.0)))
        .isRight
    )
    assert(SecondOrderCascade.fromSosMatrix(IArray(IArray(1.0))).isLeft)

  test("Fir requires causal kernel"):
    val centered = Kernel.centeredOdd(Vec(0.25, 0.5, 0.25)).orThrow
    assert(Fir(centered).isLeft)
    val fir = Fir.causal(Vec(1.0, 2.0, 1.0)).orThrow
    assertEquals(fir.stateLength, 2)

  test("FIR processInto length check"):
    val fir = Fir.causal(Vec(1.0, 0.0)).orThrow
    val runner = fir.newRunner()
    val err = runner.processInto(Vec(1.0, 2.0), gale.linalg.MutableDVec.zeros(1))
    assertEquals(err, Left(SignalError.LengthMismatch(2, 1)))

  test("FIR rejects and preserves invalid state lengths"):
    val fir = Fir.causal(Vec(1.0, 0.0, -0.5)).orThrow
    val invalid = FilterState.from(Vec(0.0))
    assertEquals(
      fir.newRunner(invalid),
      Left(SignalError.InvalidFilterState(expected = 2, actual = 1))
    )
    val runner = fir.newRunner()
    assertEquals(
      runner.restore(invalid),
      Left(SignalError.InvalidFilterState(expected = 2, actual = 1))
    )

  test("SOS processes impulse with two sections"):
    val bq = Biquad(1.0, 0.0, 0.0, 0.0, 0.0).orThrow // identity
    val sos = SecondOrderCascade.of(bq, bq).orThrow
    val x = Vec(1.0, 0.0, 0.0, 0.0)
    val y = sos.process(x).orThrow
    assertEqualsDouble(y(0), 1.0, 1e-15)
    assertEqualsDouble(y(1), 0.0, 1e-15)

  test("SOS runner supports state restore and detects invalid states"):
    val bq = Biquad(1.0, 0.5, 0.25, -0.2, 0.1).orThrow
    val cascade = SecondOrderCascade.of(bq).orThrow
    val runner = cascade.newRunner()
    val _ = runner.process(Vec(1.0, 0.0)).orThrow
    val state = runner.snapshot
    val continued = runner.process(Vec(0.0, 0.0)).orThrow
    runner.restore(state).orThrow
    assertClose(runner.process(Vec(0.0, 0.0)).orThrow, continued)
    assert(
      SecondOrderCascade(gain = Double.NaN, IArray.empty[Biquad]).isLeft,
      "non-finite cascade gain must be rejected"
    )
    assert(cascade.newRunner(FilterState.from(Vec(0.0))).isLeft)
    assert(runner.processInto(Vec(1.0), gale.linalg.MutableDVec.zeros(2)).isLeft)

  test("TF normalizes a0"):
    val tf = DigitalTransferFunction(Vec(2.0, 0.0), Vec(2.0, -1.0)).orThrow
    assertEqualsDouble(tf.a(0), 1.0, 1e-15)
    assertEqualsDouble(tf.b(0), 1.0, 1e-15)

  test("transfer-function validation, runner state, and reset are explicit"):
    assert(DigitalTransferFunction(DVec.zeros(0), Vec(1.0)).isLeft)
    assert(DigitalTransferFunction(Vec(1.0), DVec.zeros(0)).isLeft)
    assert(DigitalTransferFunction(Vec(1.0), Vec(0.0)).isLeft)
    assert(DigitalTransferFunction(Vec(Double.NaN), Vec(1.0)).isLeft)
    assert(DigitalTransferFunction(Vec(1.0), Vec(1.0, Double.PositiveInfinity)).isLeft)

    val tf = DigitalTransferFunction(Vec(1.0), Vec(1.0, -0.5)).orThrow
    val runner = tf.newRunner()
    val impulse = Vec(1.0, 0.0, 0.0, 0.0)
    val first = runner.process(impulse).orThrow
    assertEqualsDouble(first(0), 1.0, 1e-15)
    assertEqualsDouble(first(1), 0.5, 1e-15)
    val state = runner.snapshot
    val continued = runner.process(Vec(0.0, 0.0)).orThrow
    assertEquals(runner.restore(state), Right(()))
    assertClose(runner.process(Vec(0.0, 0.0)).orThrow, continued)
    runner.reset()
    assertEqualsDouble(runner.process(Vec(0.0)).orThrow(0), 0.0, 1e-15)
    assertEquals(
      tf.newRunner(FilterState.from(DVec.zeros(0))),
      Left(SignalError.InvalidFilterState(expected = 1, actual = 0))
    )
    assertEquals(
      runner.processInto(Vec(1.0), gale.linalg.MutableDVec.zeros(2)),
      Left(SignalError.LengthMismatch(expected = 1, actual = 2))
    )

  test("ZeroPhase OddPad returns same length"):
    val fir = Fir.causal(Vec(0.25, 0.5, 0.25)).orThrow
    val x = Vec.tabulate(64)(i => math.sin(0.2 * i))
    val y = ZeroPhase.filter(x, fir, EdgeTreatment.OddPad).orThrow
    assertEquals(y.length, x.length)

  test("ZeroPhase supports every implemented FIR edge treatment"):
    val fir = Fir.causal(Vec(0.25, 0.5, 0.25)).orThrow
    val x = Vec.tabulate(64)(i => math.sin(0.2 * i) + 0.01 * i)
    List(
      EdgeTreatment.OddPad,
      EdgeTreatment.EvenPad,
      EdgeTreatment.ConstantPad(-0.25),
      EdgeTreatment.Gustafsson
    ).foreach { edge =>
      val y = ZeroPhase.filter(x, fir, edge).orThrow
      assertEquals(y.length, x.length, clue = edge.toString)
      var i = 0
      while i < y.length do
        assert(y(i).isFinite, s"$edge produced a non-finite sample at $i")
        i += 1
    }

  test("ZeroPhase transfer-function, cascade, and Signal overloads preserve their contracts"):
    val fir = Fir.causal(Vec(0.25, 0.5, 0.25)).orThrow
    val tf = DigitalTransferFunction(fir.taps, Vec(1.0)).orThrow
    val x = Vec.tabulate(64)(i => math.cos(0.15 * i) + 0.02 * i)
    val edges = List(EdgeTreatment.OddPad, EdgeTreatment.EvenPad, EdgeTreatment.ConstantPad(0.5))
    edges.foreach { edge =>
      assertClose(
        ZeroPhase.filter(x, fir, edge).orThrow,
        ZeroPhase.filter(x, tf, edge).orThrow
      )
    }
    val fs = SampleRate.hertz(40.0).orThrow
    val signal = Signal(x, Sampling(fs, Seconds.of(0.25).orThrow)).orThrow
    val filtered = ZeroPhase.filter(signal, tf, EdgeTreatment.OddPad).orThrow
    assertEquals(filtered.sampleRate, fs)
    assertEqualsDouble(filtered.start.value, signal.start.value, 1e-15)

    val bq = Biquad(1.0, 0.2, 0.1, -0.15, 0.05).orThrow
    val cascade = SecondOrderCascade.of(bq).orThrow
    edges.foreach { edge =>
      val y = ZeroPhase.filter(x, cascade, edge).orThrow
      assertEquals(y.length, x.length)
      assert(y(0).isFinite, clue = edge.toString)
    }

  test("ZeroPhase reports unsupported SOS Gustafsson and too-short padded input"):
    val bq = Biquad(1.0, 0.0, 0.0, 0.0, 0.0).orThrow
    val sos = SecondOrderCascade.of(bq).orThrow
    assert(ZeroPhase.filter(Vec(1.0, 2.0, 3.0, 4.0), sos, EdgeTreatment.Gustafsson).isLeft)
    val fir = Fir.causal(Vec(0.25, 0.5, 0.25)).orThrow
    assert(ZeroPhase.filter(Vec(1.0, 2.0, 3.0), fir, EdgeTreatment.EvenPad).isLeft)

  private def assertClose(actual: DVec, expected: DVec): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assertEqualsDouble(actual(i), expected(i), 1e-15, clue = s"i=$i")
      i += 1
