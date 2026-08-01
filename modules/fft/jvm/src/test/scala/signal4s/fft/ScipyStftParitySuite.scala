package signal4s.fft

import gale.linalg.DVec
import signal4s.*
import java.nio.file.{Files, Paths}

class ScipyStftParitySuite extends munit.FunSuite:

  private val smoke = Paths.get("fixtures", "data", "smoke")

  test("ShortTimeFFT Hann periodic hop16 analyze + dual + round-trip"):
    val text = read("smoke.stft_hann_periodic_hop16")
    val signal = DVec.fromSeq(arrayAfter(text, "inputs", "signal"))
    val expRe = arrayAfter(text, "expected", "stft_real")
    val expIm = arrayAfter(text, "expected", "stft_imag")
    val expDual = DVec.fromSeq(arrayAfter(text, "expected", "dual"))
    val expRt = DVec.fromSeq(arrayAfter(text, "expected", "roundtrip"))
    val fs = SampleRate.hertz(1000.0).orThrow
    val plan = StftPlan.hannPeriodic(32, fs, hop = Some(16), nfft = Some(32)).orThrow
    assert(plan.dual.isDefined)
    assertClose(plan.dual.get.taps, expDual, 1e-9, 1e-11)

    val tf = plan.analyze(signal).orThrow
    assertEquals(tf.frameCount * tf.binCount, expRe.length)
    var i = 0
    var f = 0
    while f < tf.frameCount do
      val col = tf.columns(f)
      var k = 0
      while k < col.length do
        val tol = 1e-11 + 1e-9 * math.hypot(expRe(i), expIm(i))
        assert(
          math.hypot(col.real(k) - expRe(i), col.imaginary(k) - expIm(i)) <= tol,
          s"frame=$f bin=$k"
        )
        i += 1
        k += 1
      f += 1

    val y = plan.synthesize(tf, signal.length).orThrow
    assertClose(y, expRt, 1e-9, 1e-11)
    assertClose(y, signal, 1e-9, 1e-10)

  private def read(id: String): String =
    val path = smoke.resolve(s"$id.json")
    assume(Files.isRegularFile(path), s"missing $path")
    Files.readString(path)

  private def arrayAfter(text: String, section: String, name: String): List[Double] =
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
