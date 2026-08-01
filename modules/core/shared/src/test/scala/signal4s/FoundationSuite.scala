package signal4s

import gale.linalg.Vec

class FoundationSuite extends munit.FunSuite:

  test("SampleRate rejects non-positive and non-finite values"):
    assert(SampleRate.hertz(0.0).isLeft)
    assert(SampleRate.hertz(-1.0).isLeft)
    assert(SampleRate.hertz(Double.NaN).isLeft)
    assert(SampleRate.hertz(Double.PositiveInfinity).isLeft)
    assertEquals(SampleRate.hertz(1000.0).map(_.hertz), Right(1000.0))

  test("Frequency rejects negative and non-finite values"):
    assert(Frequency.hertz(-0.1).isLeft)
    assert(Frequency.hertz(Double.NaN).isLeft)
    assertEquals(Frequency.hertz(0.0).map(_.hertz), Right(0.0))
    assertEquals(Frequency.hertz(40.0).map(_.hertz), Right(40.0))

  test("normalized frequency conversions require sample rate and respect Nyquist"):
    val fs = SampleRate.hertz(1000.0).orThrow
    val ok = Frequency.hertz(40.0).orThrow
    val cycles = ok.cyclesPerSample(at = fs).orThrow
    val radians = ok.radiansPerSample(at = fs).orThrow
    assertEqualsDouble(cycles.value, 0.04, 1e-15)
    assertEqualsDouble(radians.value, 2.0 * math.Pi * 0.04, 1e-15)

    val above = Frequency.hertz(600.0).orThrow
    assertEquals(
      above.cyclesPerSample(at = fs),
      Left(SignalError.FrequencyAboveNyquist(above, fs))
    )

  test("Kernel.causal places origin at index 0"):
    val kernel = Kernel.causal(Vec(1.0, 2.0, 3.0)).orThrow
    assertEquals(kernel.zeroLagIndex, 0)
    assertEquals(kernel.minLag, 0)
    assertEquals(kernel.maxLag, 2)

  test("Kernel.at rejects out-of-range origins and empty taps"):
    assertEquals(Kernel.at(Vec(), 0), Left(SignalError.EmptyKernel))
    assertEquals(
      Kernel.at(Vec(1.0, 2.0), zeroLagIndex = 2),
      Left(SignalError.InvalidKernelOrigin(2, 2))
    )
    assertEquals(
      Kernel.at(Vec(1.0, 2.0), zeroLagIndex = -1),
      Left(SignalError.InvalidKernelOrigin(-1, 2))
    )

  test("Kernel.centeredOdd rejects even lengths and centers odd lengths"):
    assertEquals(
      Kernel.centeredOdd(Vec(1.0, 2.0, 3.0, 4.0)),
      Left(SignalError.KernelNotOddLength(4))
    )
    val odd = Kernel.centeredOdd(Vec(0.2, 0.6, 0.2)).orThrow
    assertEquals(odd.zeroLagIndex, 1)
    assertEquals(odd.minLag, -1)
    assertEquals(odd.maxLag, 1)

  test("Signal carries sampling and rejects empty samples"):
    val fs = SampleRate.hertz(1000.0).orThrow
    val sampling = Sampling(fs, Seconds.Zero)
    assertEquals(Signal(Vec(), sampling), Left(SignalError.EmptySignal))
    val signal = Signal(Vec.tabulate(8)(_.toDouble), sampling).orThrow
    assertEquals(signal.length, 8)
    assertEquals(signal.sampleRate.hertz, 1000.0)
    assertEqualsDouble(signal.timeAt(3).value, 0.003, 1e-15)

  test("kernel output start offset uses zero-lag origin"):
    val fs = SampleRate.hertz(100.0).orThrow
    val kernel = Kernel.at(Vec(1.0, 2.0, 3.0), zeroLagIndex = 2).orThrow
    assertEqualsDouble(kernel.outputStartOffset(fs).value, -0.02, 1e-15)
