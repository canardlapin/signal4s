package signal4s.testkit

import java.nio.file.{Files, Path, Paths}
import scala.jdk.CollectionConverters.*

/** JVM test-only reader for committed SciPy/NumPy fixture JSON.
  *
  * It validates the invariant metadata shared by every fixture before exposing
  * typed scalar and vector accessors, so parity suites cannot accidentally
  * ignore a changed pin or tolerance.
  */
final class SciPyFixture private (private val json: ujson.Obj):

  def id: String = string("id")

  def operation: String = string("operation")

  def array(section: String, name: String): Vector[Double] =
    objectAt(section)(name).arr.iterator.map(_.num).toVector

  def number(section: String, name: String): Double =
    objectAt(section)(name).num

  def int(section: String, name: String): Int =
    number(section, name).toInt

  def string(name: String): String =
    json(name).str

  def string(section: String, name: String): String =
    objectAt(section)(name).str

  def rtol: Double = number("tolerance", "rtol")

  def atol: Double = number("tolerance", "atol")

  private def objectAt(name: String): ujson.Obj =
    json(name).obj

object SciPyFixture:

  val scipyVersion = "1.15.3"
  val numpyVersion = "2.2.6"

  private val smokeRoot = Paths.get("fixtures", "data", "smoke")
  private val required = Set("id", "operation", "versions", "params", "inputs", "expected", "tolerance")

  def load(id: String): SciPyFixture =
    loadPath(smokeRoot.resolve(s"$id.json"))

  def loadPath(path: Path): SciPyFixture =
    require(Files.isRegularFile(path), s"missing fixture $path")
    val fixture = new SciPyFixture(ujson.read(Files.readString(path)).obj)
    validate(fixture, path)
    fixture

  def allSmoke(): List[SciPyFixture] =
    val paths = Files.list(smokeRoot)
    try
      paths.iterator.asScala
        .filter(path => path.getFileName.toString.endsWith(".json"))
        .toList
        .sortBy(_.getFileName.toString)
        .map(loadPath)
    finally paths.close()

  private def validate(fixture: SciPyFixture, path: Path): Unit =
    val keys = fixture.json.value.keySet
    require(required.subsetOf(keys), s"fixture $path is missing ${required.diff(keys)}")
    require(
      fixture.string("versions", "scipy") == scipyVersion,
      s"fixture $path has unexpected SciPy pin"
    )
    require(
      fixture.string("versions", "numpy") == numpyVersion,
      s"fixture $path has unexpected NumPy pin"
    )
    require(fixture.rtol >= 0.0 && fixture.atol >= 0.0, s"fixture $path has negative tolerance")
