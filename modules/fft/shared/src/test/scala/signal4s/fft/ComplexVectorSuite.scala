package signal4s.fft

import gale.linalg.{MutableDVec, Vec}
import signal4s.*

class ComplexVectorSuite extends munit.FunSuite:

  test("complex scalar arithmetic and polar form obey the usual identities"):
    val a = Complex(3.0, 4.0)
    val b = Complex(1.0, -2.0)
    assertEquals(a + b, Complex(4.0, 2.0))
    assertEquals(a - b, Complex(2.0, 6.0))
    assertEquals(a * b, Complex(11.0, -2.0))
    assertEquals(a * 0.5, Complex(1.5, 2.0))
    assertEquals(a.conjugate, Complex(3.0, -4.0))
    assertEqualsDouble(a.abs, 5.0, 1e-15)
    val polar = Complex.fromPolar(2.0, math.Pi / 2.0)
    assertEqualsDouble(polar.real, 0.0, 1e-15)
    assertEqualsDouble(polar.imaginary, 2.0, 1e-15)

  test("split-complex vectors update, scale, conjugate, and copy independently"):
    val original = ComplexVector(Vec(1.0, -2.0), Vec(3.0, 4.0)).orThrow
    assertEquals(original(1), Complex(-2.0, 4.0))
    val updated = original.updated(0, Complex(5.0, -1.0))
    assertEquals(original(0), Complex(1.0, 3.0))
    assertEquals(updated(0), Complex(5.0, -1.0))
    assertEquals(updated.scale(2.0)(1), Complex(-4.0, 8.0))
    assertEquals(updated.conjugate(0), Complex(5.0, 1.0))
    val copied = updated.copy
    assertEquals(copied(0), updated(0))
    assertEquals(copied(1), updated(1))

  test("constructors expose shape and interleaving errors"):
    assertEquals(ComplexVector.zeros(-1), Left(SignalError.InvalidInputLength(-1)))
    assertEquals(
      ComplexVector(Vec(1.0), Vec(2.0, 3.0)),
      Left(SignalError.LengthMismatch(expected = 1, actual = 2))
    )
    assert(
      ComplexVector.fromInterleaved(Vec(1.0, 2.0, 3.0)).isLeft,
      "odd interleaved storage must be rejected"
    )
    assertEquals(ComplexVector.tabulate(-1)(_ => Complex.Zero), Left(SignalError.InvalidInputLength(-1)))

  test("interleaved and tabulated constructors preserve each component"):
    val interleaved = ComplexVector.fromInterleaved(Vec(1.0, -1.0, 2.0, -2.0)).orThrow
    assertEquals(interleaved.length, 2)
    assertEquals(interleaved(0), Complex(1.0, -1.0))
    assertEquals(interleaved(1), Complex(2.0, -2.0))

    val tabulated = ComplexVector.tabulate(3)(i => Complex(i, -i)).orThrow
    assertEquals(tabulated(2), Complex(2.0, -2.0))
    val dest = MutableDVec.zeros(3)
    assertEquals(tabulated.copyRealTo(dest), Right(()))
    assertEqualsDouble(dest(2), 2.0, 1e-15)
    assertEquals(
      tabulated.copyRealTo(MutableDVec.zeros(2)),
      Left(SignalError.LengthMismatch(expected = 3, actual = 2))
    )
