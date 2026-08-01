package signal4s.design

import gale.linalg.DVec
import signal4s.*
import signal4s.fft.{WindowConvention, WindowSpec}
import java.nio.file.{Files, Paths}

class ScipyDesignParitySuite extends munit.FunSuite:

  private val smoke = Paths.get("fixtures", "data", "smoke")
  private val fs = SampleRate.hertz(1000.0).orThrow

  test("firwin lowpass Hamming 11"):
    val text = read("smoke.firwin_lowpass_hamming_11")
    val expected = DVec.fromSeq(array(text, "expected", "taps"))
    val d = FirDesign
      .lowPass(
        11,
        Frequency.hertz(100.0).orThrow,
        fs,
        WindowSpec.Hamming(11, WindowConvention.Symmetric)
      )
      .orThrow
    assertClose(d.fir.taps, expected, 1e-10, 1e-12)

  test("kaiserord"):
    val text = read("smoke.kaiserord_40db_50hz")
    val (n, beta) = FirDesign.kaiserOrder(40.0, 50.0, fs).orThrow
    assertEquals(n, intField(text, "numtaps"))
    assertEqualsDouble(beta, numberAfter(text, "beta"), 1e-12)

  test("butterworth ZPK + impulse/step/freqz"):
    val text = read("smoke.butter_lowpass4_100hz")
    val d = Butterworth.lowPass(4, Frequency.hertz(100.0).orThrow, fs).orThrow
    assertEqualsDouble(d.zpk.gain, numberAfter(text, "zpk_gain"), 1e-12)

    val expPoleRe = array(text, "expected", "zpk_poles_real").sorted
    val actPoleRe = d.zpk.poles.map(_.real).toList.sorted
    assertEquals(actPoleRe.length, expPoleRe.length)
    var i = 0
    while i < actPoleRe.length do
      assertEqualsDouble(actPoleRe(i), expPoleRe(i), 1e-9)
      i += 1

    val imp = Response.impulse(d.sos, 16).orThrow
    assertClose(imp, DVec.fromSeq(array(text, "expected", "impulse")), 1e-8, 1e-10)

    val step = Response.step(d.sos, 16).orThrow
    assertClose(step, DVec.fromSeq(array(text, "expected", "step")), 1e-8, 1e-10)

    val fr = Response.sosFreqz(d.sos, fs, 16).orThrow
    assertClose(fr.magnitude, DVec.fromSeq(array(text, "expected", "sosfreqz_mag")), 1e-7, 1e-9)

    val tf = Transfer.fromZerosPolesGain(d.zpk).orThrow
    val frTf = Response.freqz(tf, fs, 16).orThrow
    assertClose(frTf.magnitude, DVec.fromSeq(array(text, "expected", "freqz_mag")), 1e-7, 1e-9)

  private def read(id: String): String =
    val path = smoke.resolve(s"$id.json")
    assume(Files.isRegularFile(path), s"missing $path")
    Files.readString(path)

  private def array(text: String, section: String, name: String): List[Double] =
    val sec = text.indexOf(s""""$section"""")
    assert(sec >= 0)
    val nameIdx = text.indexOf(s""""$name"""", sec)
    assert(nameIdx >= 0)
    val bracket = text.indexOf('[', nameIdx)
    val end = text.indexOf(']', bracket)
    val body = text.substring(bracket + 1, end).trim
    if body.isEmpty then Nil else body.split(',').toList.map(_.trim.toDouble)

  private def numberAfter(text: String, name: String): Double =
    raw""""$name"\s*:\s*([0-9eE.+-]+)""".r.findFirstMatchIn(text).get.group(1).toDouble

  private def intField(text: String, name: String): Int =
    numberAfter(text, name).toInt

  private def assertClose(a: DVec, b: DVec, rtol: Double, atol: Double): Unit =
    assertEquals(a.length, b.length)
    var i = 0
    while i < a.length do
      val tol = atol + rtol * math.max(math.abs(a(i)), math.abs(b(i)))
      assert(math.abs(a(i) - b(i)) <= tol, s"i=$i ${a(i)} vs ${b(i)}")
      i += 1
