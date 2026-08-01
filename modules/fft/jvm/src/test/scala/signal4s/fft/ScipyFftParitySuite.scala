package signal4s.fft

import gale.linalg.DVec
import signal4s.*
import signal4s.testkit.SciPyFixture

class ScipyFftParitySuite extends munit.FunSuite:

  test("fft backward sine on bin"):
    checkComplex("smoke.fft_forward_sine_on_bin", FftNormalization.Backward)

  test("fft ortho length 16"):
    checkComplex("smoke.fft_forward_norm_ortho_16", FftNormalization.Orthonormal)

  test("fft forward normalization length 16"):
    checkComplex("smoke.fft_forward_norm_forward_16", FftNormalization.Forward)

  test("fft mixed-radix 12"):
    checkComplex("smoke.fft_forward_mixed_radix_12", FftNormalization.Backward)

  test("fft Bluestein prime 17"):
    checkComplex("smoke.fft_forward_prime_17", FftNormalization.Backward)

  test("rfft backward"):
    checkReal("smoke.rfft_forward_sine_on_bin", FftNormalization.Backward)

  test("rfft forward norm prime 17"):
    checkReal("smoke.rfft_forward_prime_17", FftNormalization.Forward)

  test("rfft ortho"):
    checkReal("smoke.rfft_forward_norm_ortho_32", FftNormalization.Orthonormal)

  test("fft inverse backward"):
    checkComplexInverse("smoke.fft_inverse_backward_32", FftNormalization.Backward)

  test("fft inverse forward"):
    checkComplexInverse("smoke.fft_inverse_forward_16", FftNormalization.Forward)

  test("rfft inverse forward"):
    checkRealInverse("smoke.rfft_inverse_forward_17", FftNormalization.Forward)

  test("rfft inverse ortho"):
    checkRealInverse("smoke.rfft_inverse_ortho_32", FftNormalization.Orthonormal)

  private def checkComplex(id: String, norm: FftNormalization): Unit =
    val fixture = SciPyFixture.load(id)
    val n = fixture.int("params", "length")
    val inRe = fixture.array("inputs", "real")
    val inIm = fixture.array("inputs", "imag")
    val expRe = fixture.array("expected", "real")
    val expIm = fixture.array("expected", "imag")
    assertEquals(inRe.length, n)
    val plan = FftPlan(n, norm).orThrow
    val x = ComplexVector(DVec.fromSeq(inRe), DVec.fromSeq(inIm)).orThrow
    val y = plan.forward(x).orThrow
    assertEquals(y.length, expRe.length)
    var i = 0
    while i < y.length do
      val tol = fixture.atol + fixture.rtol * math.max(math.abs(expRe(i)), math.abs(expIm(i)))
      assert(math.hypot(y.real(i) - expRe(i), y.imaginary(i) - expIm(i)) <= tol, s"$id@$i")
      i += 1

  private def checkReal(id: String, norm: FftNormalization): Unit =
    val fixture = SciPyFixture.load(id)
    val n = fixture.int("params", "length")
    val signal = DVec.fromSeq(fixture.array("inputs", "signal"))
    val expRe = fixture.array("expected", "real")
    val expIm = fixture.array("expected", "imag")
    val plan = RealFftPlan(n, norm).orThrow
    val y = plan.forward(signal).orThrow
    assertEquals(y.bins.length, expRe.length)
    var i = 0
    while i < y.bins.length do
      val tol = fixture.atol + fixture.rtol * math.max(math.abs(expRe(i)), math.abs(expIm(i)))
      assert(
        math.hypot(y.bins.real(i) - expRe(i), y.bins.imaginary(i) - expIm(i)) <= tol,
        s"$id@$i"
      )
      i += 1

  private def checkComplexInverse(id: String, norm: FftNormalization): Unit =
    val fixture = SciPyFixture.load(id)
    val n = fixture.int("params", "length")
    val spectrum = ComplexVector(
      DVec.fromSeq(fixture.array("inputs", "real")),
      DVec.fromSeq(fixture.array("inputs", "imag"))
    ).orThrow
    val expectedRe = fixture.array("expected", "real")
    val expectedIm = fixture.array("expected", "imag")
    val actual = FftPlan(n, norm).orThrow.inverse(spectrum).orThrow
    var i = 0
    while i < n do
      val tol = fixture.atol + fixture.rtol * math.max(math.abs(expectedRe(i)), math.abs(expectedIm(i)))
      assert(
        math.hypot(actual.real(i) - expectedRe(i), actual.imaginary(i) - expectedIm(i)) <= tol,
        s"$id@$i"
      )
      i += 1

  private def checkRealInverse(id: String, norm: FftNormalization): Unit =
    val fixture = SciPyFixture.load(id)
    val n = fixture.int("params", "length")
    val bins = ComplexVector(
      DVec.fromSeq(fixture.array("inputs", "real")),
      DVec.fromSeq(fixture.array("inputs", "imag"))
    ).orThrow
    val axis = FrequencyAxis.normalizedRealFft(n).orThrow
    val spectrum = RealSpectrum(bins, axis, n, norm).orThrow
    val expected = fixture.array("expected", "signal")
    val actual = RealFftPlan(n, norm).orThrow.inverse(spectrum).orThrow
    var i = 0
    while i < n do
      val tol = fixture.atol + fixture.rtol * math.abs(expected(i))
      assertEqualsDouble(actual(i), expected(i), tol, clue = s"$id@$i")
      i += 1
