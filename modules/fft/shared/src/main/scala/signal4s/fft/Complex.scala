package signal4s.fft

/** Scalar complex value for small sets (zeros, poles, occasional element access).
  * Large transforms use [[ComplexVector]] split storage.
  */
final case class Complex(real: Double, imaginary: Double):
  def +(that: Complex): Complex =
    Complex(real + that.real, imaginary + that.imaginary)

  def -(that: Complex): Complex =
    Complex(real - that.real, imaginary - that.imaginary)

  def *(that: Complex): Complex =
    Complex(
      real * that.real - imaginary * that.imaginary,
      real * that.imaginary + imaginary * that.real
    )

  def *(scale: Double): Complex =
    Complex(real * scale, imaginary * scale)

  def conjugate: Complex = Complex(real, -imaginary)

  def abs: Double = math.hypot(real, imaginary)

object Complex:
  val Zero: Complex = Complex(0.0, 0.0)
  val One: Complex = Complex(1.0, 0.0)

  def fromPolar(magnitude: Double, phase: Double): Complex =
    Complex(magnitude * math.cos(phase), magnitude * math.sin(phase))
