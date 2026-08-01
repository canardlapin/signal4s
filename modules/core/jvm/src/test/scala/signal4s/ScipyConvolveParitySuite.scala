package signal4s

import gale.linalg.DVec
import java.nio.file.{Files, Paths}
import scala.util.matching.Regex

/** JVM-only SciPy fixture parity for direct convolution / correlation. */
class ScipyConvolveParitySuite extends munit.FunSuite:

  private val smoke = Paths.get("fixtures", "data", "smoke")

  test("full causal odd matches SciPy convolve full"):
    val fix = load("smoke.convolve_full_causal_odd")
    val signal = DVec.fromSeq(fix.array("inputs", "signal"))
    val kernel = Kernel.causal(DVec.fromSeq(fix.array("inputs", "kernel"))).orThrow
    val expected = DVec.fromSeq(fix.array("expected", "samples"))
    val actual = Convolution(signal, kernel, OutputRegion.Full).orThrow
    assertClose(actual, expected, fix.rtol, fix.atol)

  test("full even kernel origin0 matches SciPy"):
    val fix = load("smoke.convolve_full_even_kernel_origin0")
    val signal = DVec.fromSeq(fix.array("inputs", "signal"))
    val taps = DVec.fromSeq(fix.array("inputs", "kernel"))
    val kernel = Kernel.at(taps, zeroLagIndex = 0).orThrow
    val expected = DVec.fromSeq(fix.array("expected", "samples"))
    val actual = Convolution(signal, kernel, OutputRegion.Full).orThrow
    assertClose(actual, expected, fix.rtol, fix.atol)

  test("valid causal odd matches SciPy convolve valid"):
    val fix = load("smoke.convolve_valid_causal_odd")
    val signal = DVec.fromSeq(fix.array("inputs", "signal"))
    val kernel = Kernel.causal(DVec.fromSeq(fix.array("inputs", "kernel"))).orThrow
    val expected = DVec.fromSeq(fix.array("expected", "samples"))
    val actual = Convolution(signal, kernel, OutputRegion.Valid).orThrow
    assertClose(actual, expected, fix.rtol, fix.atol)

  test("SciPy same (odd) maps to centeredOdd + Input(Zero)"):
    val fix = load("smoke.convolve_same_odd_as_centered_input_zero")
    val signal = DVec.fromSeq(fix.array("inputs", "signal"))
    val kernel = Kernel.centeredOdd(DVec.fromSeq(fix.array("inputs", "kernel"))).orThrow
    val expected = DVec.fromSeq(fix.array("expected", "samples"))
    val actual =
      Convolution(signal, kernel, OutputRegion.Input(Boundary.Zero)).orThrow
    assertClose(actual, expected, fix.rtol, fix.atol)

  test("correlation raw matches signal4s-convention fixture"):
    val fix = load("smoke.correlate_raw_impulse_odd")
    val x = DVec.fromSeq(fix.array("inputs", "x"))
    val y = DVec.fromSeq(fix.array("inputs", "y"))
    val expected = DVec.fromSeq(fix.array("expected", "values"))
    val firstLag = fix.int("expected", "first_lag")
    val actual = Correlate(x, y).orThrow
    assertEquals(actual.lags.first, firstLag)
    assertClose(actual.values, expected, fix.rtol, fix.atol)

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
      // Locate "name": [ ... ] after the section key; fixtures are shallow.
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

    def int(section: String, name: String): Int =
      val sectionIdx = text.indexOf(s""""$section"""")
      assert(sectionIdx >= 0, s"section $section not found")
      val nameIdx = text.indexOf(s""""$name"""", sectionIdx)
      assert(nameIdx >= 0, s"field $section.$name not found")
      val colon = text.indexOf(':', nameIdx)
      val slice = text.substring(colon + 1, colon + 24).trim
      val token = slice.takeWhile(c => c == '-' || c.isDigit)
      token.toInt

    private def numberField(name: String): Double =
      val pattern = raw""""$name"\s*:\s*([0-9eE.+-]+)""".r
      pattern.findFirstMatchIn(text) match
        case Some(m) => m.group(1).toDouble
        case None    => fail(s"number $name not found")
