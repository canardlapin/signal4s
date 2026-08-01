package signal4s.fft

import gale.linalg.DVec
import signal4s.*
import signal4s.filter.SecondOrderCascade
import java.nio.file.{Files, Paths}

/** End-to-end SciPy pipeline: sosfilt → Welch (E10 audit). */
class ScipyE2eParitySuite extends munit.FunSuite:

  private val smoke = Paths.get("fixtures", "data", "smoke")

  test("butter SOS → Welch density matches SciPy e2e fixture"):
    val text = read("smoke.e2e_butter_welch")
    val signal = DVec.fromSeq(array(text, "inputs", "signal"))
    val sosFlat = array(text, "inputs", "sos_flat")
    val nSec = sosFlat.length / 6
    val rows = IArray.tabulate(nSec)(s => IArray.tabulate(6)(k => sosFlat(s * 6 + k)))
    val sos = SecondOrderCascade.fromSosMatrix(rows).orThrow
    val filtered = sos.process(signal).orThrow
    assertClose(filtered, DVec.fromSeq(array(text, "expected", "filtered")), 1e-8, 1e-10)

    val fs = SampleRate.hertz(1000.0).orThrow
    val win = Window.fromSpec(WindowSpec.Hann(64, WindowConvention.Periodic)).orThrow
    val plan = WelchPlan(
      win,
      hop = 32,
      sampleRate = fs,
      nfft = 64,
      detrend = Detrend.Mean,
      scaling = SpectralScaling.Density
    ).orThrow
    val psd = plan.estimate(filtered).orThrow
    assertClose(psd.power, DVec.fromSeq(array(text, "expected", "power")), 1e-8, 1e-10)

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

  private def assertClose(a: DVec, b: DVec, rtol: Double, atol: Double): Unit =
    assertEquals(a.length, b.length)
    var i = 0
    while i < a.length do
      val tol = atol + rtol * math.max(math.abs(a(i)), math.abs(b(i)))
      assert(math.abs(a(i) - b(i)) <= tol, s"i=$i ${a(i)} vs ${b(i)}")
      i += 1
