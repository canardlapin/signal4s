package signal4s

import gale.linalg.DVec
import java.nio.file.{Files, Paths}

/** Forward operator must match SciPy convolve fixtures (same as E2). */
class ScipyOperatorParitySuite extends munit.FunSuite:

  private val smoke = Paths.get("fixtures", "data", "smoke")

  test("operator Full matches SciPy full causal odd fixture"):
    val text = Files.readString(smoke.resolve("smoke.convolve_full_causal_odd.json"))
    val signal = DVec.fromSeq(arrayAfter(text, "signal"))
    val taps = DVec.fromSeq(arrayAfter(text, "kernel"))
    val expected = DVec.fromSeq(arrayAfter(text, "samples"))
    val kernel = Kernel.causal(taps).orThrow
    val op = Convolution.operator(kernel, signal.length, OutputRegion.Full).orThrow
    val actual = op(signal)
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assertEqualsDouble(actual(i), expected(i), 1e-12)
      i += 1

  test("operator Valid matches SciPy valid fixture"):
    val text = Files.readString(smoke.resolve("smoke.convolve_valid_causal_odd.json"))
    val signal = DVec.fromSeq(arrayAfter(text, "signal"))
    val taps = DVec.fromSeq(arrayAfter(text, "kernel"))
    val expected = DVec.fromSeq(arrayAfter(text, "samples"))
    val kernel = Kernel.causal(taps).orThrow
    val op = Convolution.operator(kernel, signal.length, OutputRegion.Valid).orThrow
    val actual = op(signal)
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assertEqualsDouble(actual(i), expected(i), 1e-12)
      i += 1

  private def arrayAfter(text: String, name: String): List[Double] =
    val nameIdx = text.indexOf(s""""$name"""")
    assert(nameIdx >= 0, s"$name not found")
    val bracket = text.indexOf('[', nameIdx)
    val end = text.indexOf(']', bracket)
    val body = text.substring(bracket + 1, end).trim
    if body.isEmpty then Nil else body.split(',').toList.map(_.trim.toDouble)
