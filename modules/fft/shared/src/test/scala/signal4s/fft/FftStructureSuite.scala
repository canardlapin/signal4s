package signal4s.fft

import gale.linalg.{DVec, Vec}
import signal4s.*

class FftStructureSuite extends munit.FunSuite:

  test("FFT workspaces clear transform buffers without replacing scratch storage"):
    val scratchRe = Array(7.0, 8.0)
    val scratchIm = Array(9.0, 10.0)
    val workspace = FftWorkspace(Array(1.0, 2.0), Array(-1.0, -2.0), scratchRe, scratchIm)
    assertEquals(workspace.length, 2)
    workspace.clear()
    assertEquals(workspace.re.toList, List(0.0, 0.0))
    assertEquals(workspace.im.toList, List(0.0, 0.0))
    assertEquals(workspace.scratchRe.toList, List(7.0, 8.0))
    assertEquals(workspace.scratchIm.toList, List(9.0, 10.0))

  test("window specifications validate every constrained family"):
    assertEquals(
      WindowSpec.Rectangular(0, WindowConvention.Periodic).validate,
      Left(SignalError.InvalidInputLength(0))
    )
    assert(WindowSpec.Kaiser(4, Double.NaN, WindowConvention.Periodic).validate.isLeft)
    assert(WindowSpec.Kaiser(4, -1.0, WindowConvention.Periodic).validate.isLeft)
    List(
      WindowSpec.Hann(1, WindowConvention.Symmetric),
      WindowSpec.Hamming(1, WindowConvention.Symmetric),
      WindowSpec.Blackman(1, WindowConvention.Symmetric)
    ).foreach(spec => assert(spec.validate.isLeft, clue = spec.toString))
    assert(WindowSpec.Hann(1, WindowConvention.Periodic).validate.isRight)

  test("time and frequency axes validate lengths, finite steps, and grids"):
    val zero = Seconds.Zero
    assert(TimeAxis(zero, Seconds.unsafe(1.0), -1).isLeft)
    assert(TimeAxis(zero, Seconds.unsafe(Double.NaN), 1).isLeft)
    val time = TimeAxis.fromFrames(3, hop = 2, SampleRate.hertz(4.0).orThrow, 1.0).orThrow
    assertEqualsDouble(time.first.value, 0.25, 1e-15)
    assertEqualsDouble(time.timeAt(2).value, 1.25, 1e-15)

    val first = Frequency.hertz(0.0).orThrow
    assert(FrequencyAxis(first, Frequency.unsafe(Double.NaN), 1).isLeft)
    assert(FrequencyAxis.realFft(0, SampleRate.hertz(4.0).orThrow).isLeft)
    val frequency = FrequencyAxis.realFft(8, SampleRate.hertz(16.0).orThrow).orThrow
    assertEquals(frequency.length, 5)
    assertEqualsDouble(frequency.frequencyAt(4).hertz, 8.0, 1e-15)

  test("frame plans extract centered boundaries, validate indices, and overlap-add"):
    val window = Window.fromSpec(WindowSpec.Rectangular(4, WindowConvention.Periodic)).orThrow
    val plan = FramePlan(4, hop = 2, FrameAlignment.Centered, Boundary.Reflect).orThrow
    assertEquals(plan.firstFrameStart, -2)
    assertEquals(plan.frameCount(5), 4)
    assertClose(plan.extract(Vec(10.0, 20.0, 30.0), 0, window).orThrow, Vec(30.0, 20.0, 10.0, 20.0))
    assert(plan.extract(Vec(1.0, 2.0), -1, window).isLeft)
    val wrongWindow = Window.fromSpec(WindowSpec.Rectangular(3, WindowConvention.Periodic)).orThrow
    assert(plan.extract(Vec(1.0, 2.0, 3.0), 0, wrongWindow).isLeft)
    val dest = Array.fill(5)(0.0)
    plan.overlapAdd(dest, Vec(1.0, 2.0, 3.0, 4.0), 1)
    assertEquals(dest.toList, List(1.0, 2.0, 3.0, 4.0, 0.0))
    assert(FramePlan(0, 1, FrameAlignment.FromStart).isLeft)
    assert(FramePlan(4, 0, FrameAlignment.FromStart).isLeft)

  test("real spectra and time-frequency values reject incompatible shapes"):
    val bins = ComplexVector.zeros(3).orThrow
    val freqs = FrequencyAxis.normalizedRealFft(4).orThrow
    assertEquals(RealSpectrum(bins, freqs, 4, FftNormalization.Backward).orThrow.length, 3)
    assert(RealSpectrum(ComplexVector.zeros(2).orThrow, freqs, 4, FftNormalization.Backward).isLeft)
    assert(RealSpectrum(bins, FrequencyAxis.normalizedRealFft(2).orThrow, 4, FftNormalization.Backward).isLeft)
    assert(RealSpectrum(bins, freqs, 0, FftNormalization.Backward).isLeft)

    val times = TimeAxis(Seconds.Zero, Seconds.unsafe(1.0), 1).orThrow
    assert(TimeFrequency(Vector.empty, times, freqs, nfft = 4, onesided = true).isLeft)
    assert(TimeFrequency(Vector(bins), times, FrequencyAxis.normalizedRealFft(2).orThrow, 4, true).isLeft)
    assert(TimeFrequency(Vector(bins), times, freqs, 0, true).isLeft)
    val value = TimeFrequency(Vector(ComplexVector(Vec(3.0, 4.0, 0.0), Vec(4.0, 3.0, 0.0)).orThrow), times, freqs, 4, true).orThrow
    assertEquals(value.frameCount, 1)
    assertEquals(value.binCount, 3)
    assertEqualsDouble(value.spectrogram.power.head(0), 25.0, 1e-15)

  test("periodogram enforces the window and signal length contract"):
    val fs = SampleRate.hertz(100.0).orThrow
    val window = Window.fromSpec(WindowSpec.Hann(4, WindowConvention.Periodic)).orThrow
    assert(Periodogram(Vec(1.0, 2.0, 3.0), fs, window, nfft = 4).isLeft)

  private def assertClose(actual: DVec, expected: DVec): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assertEqualsDouble(actual(i), expected(i), 1e-15, clue = s"i=$i")
      i += 1
