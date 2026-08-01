package signal4s.design

import gale.linalg.DVec
import signal4s.*
import java.nio.file.{Files, Paths}

class ScipyResamplePolyParitySuite extends munit.FunSuite:

  private val smoke = Paths.get("fixtures", "data", "smoke")

  List((3, 2), (2, 3), (4, 2)).foreach { case (up, down) =>
    test(s"resample_poly up=$up down=$down"):
      val text = read(s"smoke.resample_poly_up${up}_down$down")
      val x = DVec.fromSeq(array(text, "inputs", "signal"))
      val expected = DVec.fromSeq(array(text, "expected", "samples"))
      val actual = ResamplePoly(x, up, down).orThrow
      assertClose(actual, expected, 1e-8, 1e-10)
  }

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
