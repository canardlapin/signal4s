package signal4s.fft

import gale.linalg.{DVec, DVecBuilder, MutableDVec}
import signal4s.SignalError

/** Split-complex vector: parallel real/imag [[DVec]] storage without boxing. */
final case class ComplexVector private (
    real: DVec,
    imaginary: DVec
):
  def length: Int = real.length

  def apply(index: Int): Complex =
    Complex(real(index), imaginary(index))

  def updated(index: Int, value: Complex): ComplexVector =
    ComplexVector.unsafe(real.updated(index, value.real), imaginary.updated(index, value.imaginary))

  def copy: ComplexVector =
    ComplexVector.unsafe(real.copy, imaginary.copy)

  def scale(alpha: Double): ComplexVector =
    ComplexVector.unsafe(
      DVec.tabulate(length)(i => alpha * real(i)),
      DVec.tabulate(length)(i => alpha * imaginary(i))
    )

  def conjugate: ComplexVector =
    ComplexVector.unsafe(real.copy, DVec.tabulate(length)(i => -imaginary(i)))

  def copyRealTo(dest: MutableDVec): Either[SignalError, Unit] =
    if dest.length != length then Left(SignalError.LengthMismatch(length, dest.length))
    else
      dest := real
      Right(())

object ComplexVector:
  def zeros(length: Int): Either[SignalError, ComplexVector] =
    if length < 0 then Left(SignalError.InvalidInputLength(length))
    else Right(unsafe(DVec.zeros(length), DVec.zeros(length)))

  def apply(real: DVec, imaginary: DVec): Either[SignalError, ComplexVector] =
    if real.length != imaginary.length then
      Left(SignalError.LengthMismatch(real.length, imaginary.length))
    else Right(unsafe(real, imaginary))

  def fromInterleaved(values: DVec): Either[SignalError, ComplexVector] =
    if values.length % 2 != 0 then
      Left(SignalError.NumericalFailure("ComplexVector", "interleaved length must be even"))
    else
      val n = values.length / 2
      val re = DVecBuilder.zeros(n)
      val im = DVecBuilder.zeros(n)
      var i = 0
      while i < n do
        re(i) = values(2 * i)
        im(i) = values(2 * i + 1)
        i += 1
      Right(unsafe(re.result(), im.result()))

  def tabulate(length: Int)(f: Int => Complex): Either[SignalError, ComplexVector] =
    if length < 0 then Left(SignalError.InvalidInputLength(length))
    else
      val re = DVecBuilder.zeros(length)
      val im = DVecBuilder.zeros(length)
      var i = 0
      while i < length do
        val z = f(i)
        re(i) = z.real
        im(i) = z.imaginary
        i += 1
      Right(unsafe(re.result(), im.result()))

  private[fft] def unsafe(real: DVec, imaginary: DVec): ComplexVector =
    new ComplexVector(real, imaginary)
