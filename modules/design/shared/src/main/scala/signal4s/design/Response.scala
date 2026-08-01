package signal4s.design

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*
import signal4s.filter.{DigitalTransferFunction, Fir, SecondOrderCascade}
import signal4s.fft.Complex

/** Frequency / time responses for designed filters. */
object Response:

  final case class FrequencyResponse(
      frequencies: DVec,
      real: DVec,
      imaginary: DVec
  ):
    def magnitude: DVec =
      DVec.tabulate(real.length)(i => math.hypot(real(i), imaginary(i)))

    def phase: DVec =
      DVec.tabulate(real.length)(i => math.atan2(imaginary(i), real(i)))

  def freqz(
      tf: DigitalTransferFunction,
      sampleRate: SampleRate,
      nfft: Int = 512
  ): Either[SignalError, FrequencyResponse] =
    if nfft <= 0 then Left(SignalError.InvalidInputLength(nfft))
    else
      val freqs = DVecBuilder.zeros(nfft)
      val re = DVecBuilder.zeros(nfft)
      val im = DVecBuilder.zeros(nfft)
      var k = 0
      while k < nfft do
        // SciPy freqz: ω = π * k / nfft for k = 0..nfft-1
        val omega = math.Pi * k.toDouble / nfft.toDouble
        freqs(k) = sampleRate.hertz * omega / (2.0 * math.Pi)
        val h = evalPolyRatio(tf.b, tf.a, omega)
        re(k) = h.real
        im(k) = h.imaginary
        k += 1
      Right(FrequencyResponse(freqs.result(), re.result(), im.result()))

  def sosFreqz(
      sos: SecondOrderCascade,
      sampleRate: SampleRate,
      nfft: Int = 512
  ): Either[SignalError, FrequencyResponse] =
    if nfft <= 0 then Left(SignalError.InvalidInputLength(nfft))
    else
      val freqs = DVecBuilder.zeros(nfft)
      val re = DVecBuilder.zeros(nfft)
      val im = DVecBuilder.zeros(nfft)
      var k = 0
      while k < nfft do
        val omega = math.Pi * k.toDouble / nfft.toDouble
        freqs(k) = sampleRate.hertz * omega / (2.0 * math.Pi)
        var h = Complex(sos.gain, 0.0)
        var s = 0
        while s < sos.sections.length do
          val sec = sos.sections(s)
          val num = evalBiquadNum(sec.b0, sec.b1, sec.b2, omega)
          val den = evalBiquadNum(1.0, sec.a1, sec.a2, omega)
          h = h * div(num, den)
          s += 1
        re(k) = h.real
        im(k) = h.imaginary
        k += 1
      Right(FrequencyResponse(freqs.result(), re.result(), im.result()))

  def impulse(
      fir: Fir,
      n: Int
  ): Either[SignalError, DVec] =
    if n <= 0 then Left(SignalError.InvalidInputLength(n))
    else
      val x = DVec.tabulate(n)(i => if i == 0 then 1.0 else 0.0)
      fir.process(x)

  def impulse(
      sos: SecondOrderCascade,
      n: Int
  ): Either[SignalError, DVec] =
    if n <= 0 then Left(SignalError.InvalidInputLength(n))
    else
      val x = DVec.tabulate(n)(i => if i == 0 then 1.0 else 0.0)
      sos.process(x)

  def step(
      sos: SecondOrderCascade,
      n: Int
  ): Either[SignalError, DVec] =
    if n <= 0 then Left(SignalError.InvalidInputLength(n))
    else
      val x = DVec.tabulate(n)(_ => 1.0)
      sos.process(x)

  /** Approximate group delay via -dθ/dω on uniform grid (unwrap first). */
  def groupDelay(
      sos: SecondOrderCascade,
      sampleRate: SampleRate,
      nfft: Int = 512
  ): Either[SignalError, DVec] =
    sosFreqz(sos, sampleRate, nfft).map { fr =>
      val phase = unwrap(fr.phase)
      val gd = DVecBuilder.zeros(nfft)
      val dOmega = math.Pi / nfft.toDouble
      var k = 0
      while k < nfft do
        val left = if k == 0 then phase(0) else phase(k - 1)
        val right = if k == nfft - 1 then phase(k) else phase(k + 1)
        val denom = if k == 0 || k == nfft - 1 then dOmega else 2.0 * dOmega
        gd(k) = -(right - left) / denom
        k += 1
      gd.result()
    }

  private def evalPolyRatio(b: DVec, a: DVec, omega: Double): Complex =
    div(evalPoly(b, omega), evalPoly(a, omega))

  private def evalPoly(c: DVec, omega: Double): Complex =
    var re = 0.0
    var im = 0.0
    var i = 0
    while i < c.length do
      val ang = -omega * i
      re += c(i) * math.cos(ang)
      im += c(i) * math.sin(ang)
      i += 1
    Complex(re, im)

  private def evalBiquadNum(c0: Double, c1: Double, c2: Double, omega: Double): Complex =
    val e1 = Complex(math.cos(-omega), math.sin(-omega))
    val e2 = Complex(math.cos(-2 * omega), math.sin(-2 * omega))
    Complex(c0, 0.0) + e1 * c1 + e2 * c2

  private def div(a: Complex, b: Complex): Complex =
    val n = b.real * b.real + b.imaginary * b.imaginary
    Complex(
      (a.real * b.real + a.imaginary * b.imaginary) / n,
      (a.imaginary * b.real - a.real * b.imaginary) / n
    )

  private def unwrap(phase: DVec): DVec =
    val out = DVecBuilder.zeros(phase.length)
    if phase.length == 0 then out.result()
    else
      out(0) = phase(0)
      var i = 1
      while i < phase.length do
        var dp = phase(i) - phase(i - 1)
        while dp > math.Pi do dp -= 2.0 * math.Pi
        while dp < -math.Pi do dp += 2.0 * math.Pi
        out(i) = out(i - 1) + dp
        i += 1
      out.result()
