package signal4s.fft

import gale.linalg.DVec
import signal4s.*
import java.nio.file.{Files, Paths}

/** SciPy convolve fixtures via FFT and OLA paths (same cases as E2). */
class ScipyFftConvolveParitySuite extends munit.FunSuite:

  FftBackend.ensureInstalled()

  private val smoke = Paths.get("fixtures", "data", "smoke")

  test("FFT full causal odd matches SciPy"):
    check("smoke.convolve_full_causal_odd", OutputRegion.Full, ConvolutionMethod.Fft)

  test("OLA full causal odd matches SciPy"):
    check(
      "smoke.convolve_full_causal_odd",
      OutputRegion.Full,
      ConvolutionMethod.OverlapAdd(4)
    )

  test("FFT valid causal odd matches SciPy"):
    check("smoke.convolve_valid_causal_odd", OutputRegion.Valid, ConvolutionMethod.Fft)

  test("FFT same/centered Input(Zero) matches SciPy"):
    val fix = load("smoke.convolve_same_odd_as_centered_input_zero")
    val signal = DVec.fromSeq(fix.array("inputs", "signal"))
    val kernel = Kernel.centeredOdd(DVec.fromSeq(fix.array("inputs", "kernel"))).orThrow
    val expected = DVec.fromSeq(fix.array("expected", "samples"))
    val actual =
      Convolution(signal, kernel, OutputRegion.Input(Boundary.Zero), ConvolutionMethod.Fft).orThrow
    assertClose(actual, expected, fix.rtol, fix.atol)

  test("FFT full even kernel origin0 matches SciPy"):
    val fix = load("smoke.convolve_full_even_kernel_origin0")
    val signal = DVec.fromSeq(fix.array("inputs", "signal"))
    val taps = DVec.fromSeq(fix.array("inputs", "kernel"))
    val kernel = Kernel.at(taps, zeroLagIndex = 0).orThrow
    val expected = DVec.fromSeq(fix.array("expected", "samples"))
    val actual =
      Convolution(signal, kernel, OutputRegion.Full, ConvolutionMethod.Fft).orThrow
    assertClose(actual, expected, fix.rtol, fix.atol)

  private def check(id: String, region: OutputRegion, method: ConvolutionMethod): Unit =
    val fix = load(id)
    val signal = DVec.fromSeq(fix.array("inputs", "signal"))
    val kernel = Kernel.causal(DVec.fromSeq(fix.array("inputs", "kernel"))).orThrow
    val expected = DVec.fromSeq(fix.array("expected", "samples"))
    val actual = Convolution(signal, kernel, region, method).orThrow
    assertClose(actual, expected, fix.rtol, fix.atol)

  private def load(id: String): Fixture =
    val path = smoke.resolve(s"$id.json")
    assume(Files.isRegularFile(path), s"missing $path")
    Fixture(Files.readString(path))

  private def assertClose(
      actual: DVec,
      expected: DVec,
      rtol: Double,
      atol: Double
  ): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      val tol = atol + rtol * math.max(math.abs(actual(i)), math.abs(expected(i)))
      assert(
        math.abs(actual(i) - expected(i)) <= tol,
        s"index $i: ${actual(i)} vs ${expected(i)} (tol=$tol)"
      )
      i += 1

  private final class Fixture(text: String):
    val rtol: Double = numberField("rtol")
    val atol: Double = numberField("atol")

    def array(section: String, name: String): List[Double] =
      val sectionIdx = text.indexOf(s""""$section"""")
      assert(sectionIdx >= 0, s"section $section not found")
      val nameIdx = text.indexOf(s""""$name"""", sectionIdx)
      assert(nameIdx >= 0, s"field $section.$name not found")
      val bracket = text.indexOf('[', nameIdx)
      val end = text.indexOf(']', bracket)
      assert(bracket >= 0 && end > bracket, s"array $section.$name not found")
      val body = text.substring(bracket + 1, end).trim
      if body.isEmpty then Nil
      else body.split(',').toList.map(_.trim.toDouble)

    private def numberField(name: String): Double =
      val pattern = raw""""$name"\s*:\s*([0-9eE.+-]+)""".r
      pattern.findFirstMatchIn(text) match
        case Some(m) => m.group(1).toDouble
        case None    => fail(s"number $name not found")
