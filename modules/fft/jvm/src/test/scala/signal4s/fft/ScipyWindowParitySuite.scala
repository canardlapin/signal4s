package signal4s.fft

import gale.linalg.DVec
import signal4s.*
import java.nio.file.{Files, Paths}

class ScipyWindowParitySuite extends munit.FunSuite:

  private val smoke = Paths.get("fixtures", "data", "smoke")

  private val cases = List(
    ("smoke.window_hann_periodic_16", WindowSpec.Hann(16, WindowConvention.Periodic)),
    ("smoke.window_hann_symmetric_16", WindowSpec.Hann(16, WindowConvention.Symmetric)),
    ("smoke.window_hamming_periodic_16", WindowSpec.Hamming(16, WindowConvention.Periodic)),
    ("smoke.window_hamming_symmetric_16", WindowSpec.Hamming(16, WindowConvention.Symmetric)),
    ("smoke.window_blackman_periodic_16", WindowSpec.Blackman(16, WindowConvention.Periodic)),
    ("smoke.window_blackman_symmetric_16", WindowSpec.Blackman(16, WindowConvention.Symmetric)),
    ("smoke.window_kaiser_symmetric_16", WindowSpec.Kaiser(16, 5.0, WindowConvention.Symmetric)),
    ("smoke.window_kaiser_periodic_16", WindowSpec.Kaiser(16, 5.0, WindowConvention.Periodic))
  )

  cases.foreach { case (id, spec) =>
    test(s"$id matches SciPy"):
      val fix = load(id)
      val expected = DVec.fromSeq(fix.array("expected", "taps"))
      val actual = Window.fromSpec(spec).orThrow.taps
      assertClose(actual, expected, fix.rtol, fix.atol)
  }

  private def load(id: String): Fixture =
    val path = smoke.resolve(s"$id.json")
    assume(Files.isRegularFile(path), s"missing $path")
    Fixture(Files.readString(path))

  private def assertClose(actual: DVec, expected: DVec, rtol: Double, atol: Double): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      val tol = atol + rtol * math.max(math.abs(actual(i)), math.abs(expected(i)))
      assert(math.abs(actual(i) - expected(i)) <= tol, s"i=$i ${actual(i)} vs ${expected(i)}")
      i += 1

  private final class Fixture(text: String):
    val rtol: Double = numberField("rtol")
    val atol: Double = numberField("atol")
    def array(section: String, name: String): List[Double] =
      val sectionIdx = text.indexOf(s""""$section"""")
      assert(sectionIdx >= 0)
      val nameIdx = text.indexOf(s""""$name"""", sectionIdx)
      assert(nameIdx >= 0)
      val bracket = text.indexOf('[', nameIdx)
      val end = text.indexOf(']', bracket)
      val body = text.substring(bracket + 1, end).trim
      if body.isEmpty then Nil else body.split(',').toList.map(_.trim.toDouble)
    private def numberField(name: String): Double =
      raw""""$name"\s*:\s*([0-9eE.+-]+)""".r.findFirstMatchIn(text).get.group(1).toDouble
