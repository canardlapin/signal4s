package signal4s.design

import signal4s.SignalError
import signal4s.fft.Complex

/** Analog prototype transforms (SciPy `buttap`, `lp2lp_zpk`, `lp2hp_zpk`, `lp2bp_zpk`). */
private[design] object AnalogPrototype:

  /** Unit-cutoff Butterworth analog prototype (`scipy.signal.buttap`). */
  def butterworth(order: Int): Either[SignalError, ZerosPolesGain] =
    if order <= 0 then Left(SignalError.InvalidInputLength(order))
    else
      // m = -N+1, -N+3, ..., N-1; p = -exp(1j * pi * m / (2N))
      val poles = IArray.tabulate(order) { idx =>
        val m = -order + 1 + 2 * idx
        val theta = math.Pi * m.toDouble / (2.0 * order)
        Complex(-math.cos(theta), -math.sin(theta))
      }
      ZerosPolesGain(IArray.empty[Complex], poles, 1.0)

  def lp2lp(proto: ZerosPolesGain, wo: Double): Either[SignalError, ZerosPolesGain] =
    if !(wo > 0.0 && wo.isFinite) then
      Left(SignalError.NumericalFailure("lp2lp", s"wo must be positive, got $wo"))
    else
      val degree = proto.poles.length - proto.zeros.length
      val z = proto.zeros.map(r => r * wo)
      val p = proto.poles.map(r => r * wo)
      val k = proto.gain * math.pow(wo, degree.toDouble)
      ZerosPolesGain(z, p, k)

  def lp2hp(proto: ZerosPolesGain, wo: Double): Either[SignalError, ZerosPolesGain] =
    if !(wo > 0.0 && wo.isFinite) then
      Left(SignalError.NumericalFailure("lp2hp", s"wo must be positive, got $wo"))
    else
      val degree = proto.poles.length - proto.zeros.length
      val zMapped = proto.zeros.map(s => invScale(s, wo))
      val z =
        if degree > 0 then zMapped ++ IArray.tabulate(degree)(_ => Complex.Zero)
        else zMapped
      val p = proto.poles.map(s => invScale(s, wo))
      val k = proto.gain * (realProd(negate(proto.zeros)) / realProd(negate(proto.poles)))
      ZerosPolesGain(z, p, k)

  def lp2bp(proto: ZerosPolesGain, wo: Double, bw: Double): Either[SignalError, ZerosPolesGain] =
    if !(wo > 0.0 && bw > 0.0 && wo.isFinite && bw.isFinite) then
      Left(SignalError.NumericalFailure("lp2bp", "wo and bw must be positive finite"))
    else
      val degree = proto.poles.length - proto.zeros.length
      val zBp = IArray.newBuilder[Complex]
      var i = 0
      while i < proto.zeros.length do
        appendQuadratic(zBp, proto.zeros(i), wo, bw)
        i += 1
      var d = 0
      while d < degree do
        zBp += Complex.Zero
        d += 1
      val pBp = IArray.newBuilder[Complex]
      i = 0
      while i < proto.poles.length do
        appendQuadratic(pBp, proto.poles(i), wo, bw)
        i += 1
      val k = proto.gain * math.pow(bw, degree.toDouble)
      ZerosPolesGain(zBp.result(), pBp.result(), k)

  private def appendQuadratic(
      out: scala.collection.mutable.Builder[Complex, IArray[Complex]],
      root: Complex,
      wo: Double,
      bw: Double
  ): Unit =
    val half = root * (0.5 * bw)
    val disc = half * half - Complex(wo * wo, 0.0)
    val sqrtDisc = complexSqrt(disc)
    out += half + sqrtDisc
    out += half - sqrtDisc

  private def invScale(s: Complex, wo: Double): Complex =
    val n = s.real * s.real + s.imaginary * s.imaginary
    Complex(wo * s.real / n, -wo * s.imaginary / n)

  private def negate(xs: IArray[Complex]): IArray[Complex] =
    xs.map(c => Complex(-c.real, -c.imaginary))

  private def realProd(xs: IArray[Complex]): Double =
    var acc = Complex.One
    var i = 0
    while i < xs.length do
      acc = acc * xs(i)
      i += 1
    acc.real

  private def complexSqrt(z: Complex): Complex =
    val r = z.abs
    if r == 0.0 then Complex.Zero
    else
      val re = math.sqrt(math.max(0.0, (r + z.real) / 2.0))
      val im = math.copySign(math.sqrt(math.max(0.0, (r - z.real) / 2.0)), z.imaginary)
      Complex(re, im)
