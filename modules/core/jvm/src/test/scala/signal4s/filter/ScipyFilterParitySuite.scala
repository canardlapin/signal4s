package signal4s.filter

import gale.linalg.DVec
import signal4s.*
import java.nio.file.{Files, Paths}

class ScipyFilterParitySuite extends munit.FunSuite:

  private val smoke = Paths.get("fixtures", "data", "smoke")

  test("FIR runner matches SciPy lfilter"):
    val text = read("smoke.lfilter_fir_causal")
    val signal = DVec.fromSeq(arrayNamed(text, "signal"))
    val b = DVec.fromSeq(arrayNamed(text, "b"))
    val expected = DVec.fromSeq(arrayNamed(text, "samples"))
    val fir = Fir.causal(b).orThrow
    val actual = fir.process(signal).orThrow
    assertClose(actual, expected, 1e-12, 1e-14)

  test("FIR zi/zf round-trip matches SciPy"):
    val text = read("smoke.lfilter_fir_zi_zf")
    val signal = DVec.fromSeq(arrayNamed(text, "signal"))
    val b = DVec.fromSeq(arrayNamed(text, "b"))
    val zi = DVec.fromSeq(arrayNamed(text, "zi"))
    val expected = DVec.fromSeq(arrayNamed(text, "samples"))
    val zf = DVec.fromSeq(arrayNamed(text, "zf"))
    val fir = Fir.causal(b).orThrow
    val runner = fir.newRunner(FilterState.from(zi)).orThrow
    val actual = runner.process(signal).orThrow
    assertClose(actual, expected, 1e-12, 1e-14)
    assertClose(runner.snapshot.values, zf, 1e-12, 1e-14)

  test("SOS matches SciPy sosfilt butter4"):
    val text = read("smoke.sosfilt_butter4")
    val signal = DVec.fromSeq(arrayNamed(text, "signal"))
    val expected = DVec.fromSeq(arrayNamed(text, "samples"))
    val nSections = intNamed(text, "n_sections")
    val flat = arrayNamed(text, "sos_flat")
    val sos = cascadeFromFlat(flat, nSections)
    val actual = sos.process(signal).orThrow
    assertClose(actual, expected, 1e-10, 1e-12)

  test("SOS near unit circle matches SciPy"):
    val text = read("smoke.sosfilt_near_unit_circle")
    val signal = DVec.fromSeq(arrayNamed(text, "signal"))
    val expected = DVec.fromSeq(arrayNamed(text, "samples"))
    val nSections = intNamed(text, "n_sections")
    val flat = arrayNamed(text, "sos_flat")
    val sos = cascadeFromFlat(flat, nSections)
    val actual = sos.process(signal).orThrow
    assertClose(actual, expected, 1e-9, 1e-11)

  test("ZeroPhase OddPad matches SciPy filtfilt padtype=odd"):
    val text = read("smoke.filtfilt_fir_odd")
    val signal = DVec.fromSeq(arrayNamed(text, "signal"))
    val b = DVec.fromSeq(arrayNamed(text, "b"))
    val expected = DVec.fromSeq(arrayNamed(text, "samples"))
    val fir = Fir.causal(b).orThrow
    val actual = ZeroPhase.filter(signal, fir, EdgeTreatment.OddPad).orThrow
    assertClose(actual, expected, 1e-9, 1e-11)

  private def cascadeFromFlat(flat: List[Double], nSections: Int): SecondOrderCascade =
    assertEquals(flat.length, nSections * 6)
    val rows = IArray.tabulate(nSections) { s =>
      IArray.tabulate(6)(j => flat(s * 6 + j))
    }
    SecondOrderCascade.fromSosMatrix(rows).orThrow

  private def read(id: String): String =
    val path = smoke.resolve(s"$id.json")
    assume(Files.isRegularFile(path), s"missing $path")
    Files.readString(path)

  private def arrayNamed(text: String, name: String): List[Double] =
    var from = 0
    var last: Option[List[Double]] = None
    while from < text.length do
      val nameIdx = text.indexOf(s""""$name"""", from)
      if nameIdx < 0 then from = text.length
      else
        val bracket = text.indexOf('[', nameIdx)
        val end = text.indexOf(']', bracket)
        if bracket >= 0 && end > bracket then
          val body = text.substring(bracket + 1, end).trim
          last =
            Some(
              if body.isEmpty then Nil
              else body.split(',').toList.map(_.trim.toDouble)
            )
        from = nameIdx + 1
    last.getOrElse(fail(s"array $name not found"))

  private def intNamed(text: String, name: String): Int =
    val nameIdx = text.indexOf(s""""$name"""")
    assert(nameIdx >= 0, s"$name not found")
    val colon = text.indexOf(':', nameIdx)
    text.substring(colon + 1, colon + 24).trim.takeWhile(c => c == '-' || c.isDigit).toInt

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
