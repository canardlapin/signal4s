package signal4s

import gale.linalg.Vec

class ConvolutionSuite extends munit.FunSuite:

  private val x = Vec(1.0, 2.0, 3.0, 4.0)
  private val h = Kernel.causal(Vec(0.5, 1.0, 0.25)).orThrow

  test("full length and known small case"):
    val y = Convolution(x, h, OutputRegion.Full).orThrow
    assertEquals(y.length, 6)
    // hand: causal convolve
    assertEqualsDouble(y(0), 0.5, 1e-15)
    assertEqualsDouble(y(1), 2.0, 1e-15)
    assertEqualsDouble(y(2), 3.75, 1e-15)
    assertEqualsDouble(y(3), 5.5, 1e-15)
    assertEqualsDouble(y(4), 4.75, 1e-15)
    assertEqualsDouble(y(5), 1.0, 1e-15)

  test("valid region requires full overlap"):
    val y = Convolution(x, h, OutputRegion.Valid).orThrow
    assertEquals(y.length, 2)
    assertEqualsDouble(y(0), 3.75, 1e-15)
    assertEqualsDouble(y(1), 5.5, 1e-15)

  test("large causal valid output agrees with the full-convolution interior"):
    val signal = Vec.tabulate(1024)(i => math.sin(0.03 * i) + 0.01 * (i % 7))
    val kernel = Kernel.causal(Vec.tabulate(48)(i => 1.0 / (i + 1.0))).orThrow
    val actual = Convolution(signal, kernel, OutputRegion.Valid, ConvolutionMethod.Direct).orThrow
    val full = Convolution(signal, kernel, OutputRegion.Full, ConvolutionMethod.Direct).orThrow
    var i = 0
    while i < actual.length do
      assertEqualsDouble(actual(i), full(kernel.length - 1 + i), 1e-12, clue = s"valid full interior@$i")
      i += 1

  test("input-aligned zero boundary matches truncated full with origin 0"):
    val y = Convolution(x, h, OutputRegion.Input(Boundary.Zero)).orThrow
    assertEquals(y.length, x.length)
    assertEqualsDouble(y(0), 0.5, 1e-15)
    assertEqualsDouble(y(1), 2.0, 1e-15)
    assertEqualsDouble(y(2), 3.75, 1e-15)
    assertEqualsDouble(y(3), 5.5, 1e-15)

  test("large causal input agrees with the full-convolution prefix"):
    val signal = Vec.tabulate(1024)(i => math.sin(0.03 * i) + 0.01 * (i % 7))
    val kernel = Kernel.causal(Vec.tabulate(48)(i => 1.0 / (i + 1.0))).orThrow
    val actual =
      Convolution(signal, kernel, OutputRegion.Input(Boundary.Zero), ConvolutionMethod.Direct).orThrow
    val full = Convolution(signal, kernel, OutputRegion.Full, ConvolutionMethod.Direct).orThrow
    var i = 0
    while i < signal.length do
      assertEqualsDouble(actual(i), full(i), 1e-12, clue = s"causal full prefix@$i")
      i += 1

  test("nonzero kernel origin shifts Signal time axis, not full array length"):
    val fs = SampleRate.hertz(100.0).orThrow
    val signal = Signal(x, Sampling(fs, Seconds.Zero)).orThrow
    val centered = Kernel.centeredOdd(Vec(0.25, 0.5, 0.25)).orThrow
    val y = Convolution(signal, centered, OutputRegion.Full).orThrow
    assertEquals(y.length, x.length + centered.length - 1)
    assertEqualsDouble(y.start.value, -0.01, 1e-15)

  test("reflect boundary is defined at edges"):
    val k = Kernel.centeredOdd(Vec(1.0, 0.0, 0.0)).orThrow
    // y[n] = x[n - (-1)] = x[n+1] with reflect for n near end
    val y = Convolution(Vec(10.0, 20.0, 30.0), k, OutputRegion.Input(Boundary.Reflect)).orThrow
    assertEquals(y.length, 3)
    assertEqualsDouble(y(0), 20.0, 1e-15)
    assertEqualsDouble(y(1), 30.0, 1e-15)
    assertEqualsDouble(y(2), 20.0, 1e-15) // reflect: index 3 -> 1

  test("input boundaries distinguish clamp, reflect, symmetric, and constant extension"):
    val samples = Vec(10.0, 20.0, 30.0)
    val k = Kernel.centeredOdd(Vec(1.0, 2.0, 3.0)).orThrow
    val expected = List(
      Boundary.Clamp -> Vec(70.0, 100.0, 150.0),
      Boundary.Reflect -> Vec(100.0, 100.0, 140.0),
      Boundary.Symmetric -> Vec(70.0, 100.0, 150.0),
      Boundary.Constant(7.0) -> Vec(61.0, 100.0, 127.0)
    )
    expected.foreach { case (boundary, values) =>
      val actual = Convolution(samples, k, OutputRegion.Input(boundary)).orThrow
      assertClose(actual, values, clue = boundary.toString)
    }

  test("input zero honors an arbitrary kernel origin"):
    val samples = Vec(2.0, -1.0, 3.0, 4.0)
    val k = Kernel.at(Vec(1.0, -2.0, 0.5, 3.0), zeroLagIndex = 2).orThrow
    val actual = Convolution(samples, k, OutputRegion.Input(Boundary.Zero)).orThrow
    val expected = Vec.tabulate(samples.length) { i =>
      var sum = 0.0
      var j = 0
      while j < k.length do
        val source = i - j + k.zeroLagIndex
        if source >= 0 && source < samples.length then sum += k.taps(j) * samples(source)
        j += 1
      sum
    }
    assertClose(actual, expected, clue = "arbitrary-origin zero boundary")

  test("core direct plans preserve region semantics and validate input length"):
    val centered = Kernel.centeredOdd(Vec(0.25, 0.5, 0.25)).orThrow
    List(OutputRegion.Full, OutputRegion.Valid, OutputRegion.Input(Boundary.Zero)).foreach { region =>
      val plan = Convolution.plan(centered, x.length, region, ConvolutionMethod.Direct).orThrow
      assertEquals(plan.selectedMethod, ConvolutionMethod.Direct)
      assertClose(
        plan(x).orThrow,
        Convolution(x, centered, region, ConvolutionMethod.Direct).orThrow,
        clue = region.toString
      )
      assertEquals(plan(Vec(1.0)), Left(SignalError.LengthMismatch(x.length, 1)))
    }

  test("core-only convolution rejects FFT methods and invalid circular or operator shapes"):
    assert(Convolution(x, h, OutputRegion.Full, ConvolutionMethod.Fft).isLeft)
    assert(Convolution.plan(h, x.length, OutputRegion.Full, ConvolutionMethod.OverlapAdd(4)).isLeft)
    assertEquals(
      Convolution.circular(x, h, period = 0),
      Left(SignalError.InvalidPeriod(0))
    )
    val longer = Kernel.causal(Vec(1.0, 2.0, 3.0)).orThrow
    assertEquals(
      Convolution.operator(longer, inputLength = 1, OutputRegion.Valid),
      Left(SignalError.EmptyConvolution)
    )

  test("signal convolution start follows the requested output region"):
    val fs = SampleRate.hertz(100.0).orThrow
    val signal = Signal(x, Sampling(fs, Seconds.of(1.0).orThrow)).orThrow
    val centered = Kernel.centeredOdd(Vec(0.25, 0.5, 0.25)).orThrow
    val valid = Convolution(signal, centered, OutputRegion.Valid, ConvolutionMethod.Direct).orThrow
    val input = Convolution(signal, centered, OutputRegion.Input(Boundary.Zero), ConvolutionMethod.Direct).orThrow
    assertEqualsDouble(valid.start.value, 1.01, 1e-15)
    assertEqualsDouble(input.start.value, 1.0, 1e-15)

  test("circular convolution is separate and period-checked"):
    val k = Kernel.causal(Vec(1.0, 2.0)).orThrow
    val err = Convolution.circular(Vec(1.0, 2.0, 3.0), k, period = 2)
    assertEquals(err, Left(SignalError.LengthMismatch(2, 3)))
    val y = Convolution.circular(Vec(1.0, 0.0, 0.0, 0.0), k, period = 4).orThrow
    assertEqualsDouble(y(0), 1.0, 1e-15)
    assertEqualsDouble(y(1), 2.0, 1e-15)
    assertEqualsDouble(y(2), 0.0, 1e-15)

  test("correlation returns lag axis under documented convention"):
    val a = Vec(1.0, 2.0, 3.0)
    val b = Vec(0.0, 1.0, 0.5)
    val corr = Correlate(a, b).orThrow
    assertEquals(corr.lags.first, 1 - b.length)
    assertEquals(corr.lags.length, a.length + b.length - 1)
    // ℓ=0: sum_n a[n] b[n] = 1*0 + 2*1 + 3*0.5 = 3.5
    val zeroIdx = -corr.lags.first
    assertEquals(corr.lags.lagAt(zeroIdx), 0)
    assertEqualsDouble(corr.values(zeroIdx), 3.5, 1e-15)

  test("correlation normalization modes and zero-energy coefficient are explicit"):
    val a = Vec(1.0, 2.0)
    val b = Vec(3.0, 4.0)
    assertClose(
      Correlate(a, b, CorrelationNormalization.Biased).orThrow.values,
      Vec(3.0, 5.5, 2.0),
      clue = "biased"
    )
    assertClose(
      Correlate(a, b, CorrelationNormalization.Unbiased).orThrow.values,
      Vec(6.0, 5.5, 4.0),
      clue = "unbiased"
    )
    val coeff = Correlate(a, b, CorrelationNormalization.Coefficient).orThrow
    assertEqualsDouble(coeff.values(1), 11.0 / math.sqrt(125.0), 1e-15)
    val auto = Correlate.auto(a, CorrelationNormalization.Coefficient).orThrow
    assertEqualsDouble(auto.values(1), 1.0, 1e-15)
    assertClose(
      Correlate(Vec(0.0, 0.0), b, CorrelationNormalization.Coefficient).orThrow.values,
      Vec(0.0, 0.0, 0.0),
      clue = "zero-energy coefficient"
    )

  private def assertClose(actual: gale.linalg.DVec, expected: gale.linalg.DVec, clue: String): Unit =
    assertEquals(actual.length, expected.length, clue)
    var i = 0
    while i < actual.length do
      assertEqualsDouble(actual(i), expected(i), 1e-15, clue = s"$clue@$i")
      i += 1
