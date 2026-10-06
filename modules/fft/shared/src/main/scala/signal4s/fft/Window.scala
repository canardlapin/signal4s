package signal4s.fft

import gale.linalg.{DVec, DVecBuilder}
import signal4s.SignalError

/** Concrete window taps plus spectral gains.
  *
  * @param coherentGain mean(w) — amplitude / coherent processing gain
  * @param powerGain mean(w²) — power processing gain
  * @param enbwBins equivalent noise bandwidth in FFT bins:
  *                 \(N \sum w^2 / (\sum w)^2\)
  */
final case class Window private (
    taps: DVec,
    convention: WindowConvention,
    coherentGain: Double,
    powerGain: Double,
    enbwBins: Double
):
  def length: Int = taps.length

object Window:
  def fromSpec(spec: WindowSpec): Either[SignalError, Window] =
    spec.validate.flatMap { _ =>
      val taps = materialize(spec)
      fromTaps(taps, spec.convention)
    }

  def fromTaps(taps: DVec, convention: WindowConvention): Either[SignalError, Window] =
    if taps.length <= 0 then Left(SignalError.InvalidInputLength(taps.length))
    else
      val n = taps.length.toDouble
      var sum = 0.0
      var sumSq = 0.0
      var i = 0
      while i < taps.length do
        val w = taps(i)
        if !w.isFinite then
          return Left(SignalError.NumericalFailure("Window", s"nonfinite tap at index $i"))
        sum += w
        sumSq += w * w
        i += 1
      if !sum.isFinite || !sumSq.isFinite then
        return Left(SignalError.NumericalFailure("Window", "window gain exceeds finite Double capacity"))
      val coherent = sum / n
      val power = sumSq / n
      val enbw =
        if sum == 0.0 then Double.PositiveInfinity
        else n * (sumSq / sum) / sum
      Right(new Window(taps, convention, coherent, power, enbw))

  private def materialize(spec: WindowSpec): DVec =
    if spec.length == 1 then DVec.tabulate(1)(_ => 1.0)
    else spec match
      case WindowSpec.Rectangular(n, _) =>
        DVec.tabulate(n)(_ => 1.0)
      case WindowSpec.Hann(n, conv) =>
        cosineWindow(n, conv, a0 = 0.5, a1 = 0.5, a2 = 0.0)
      case WindowSpec.Hamming(n, conv) =>
        cosineWindow(n, conv, a0 = 0.54, a1 = 0.46, a2 = 0.0)
      case WindowSpec.Blackman(n, conv) =>
        cosineWindow(n, conv, a0 = 0.42, a1 = 0.5, a2 = 0.08)
      case WindowSpec.Kaiser(n, beta, conv) =>
        kaiser(n, beta, conv)

  /** SciPy cosine-sum windows (Hann/Hamming/Blackman). */
  private def cosineWindow(
      n: Int,
      convention: WindowConvention,
      a0: Double,
      a1: Double,
      a2: Double
  ): DVec =
    val out = DVecBuilder.zeros(n)
    val denom =
      convention match
        case WindowConvention.Periodic  => n.toDouble
        case WindowConvention.Symmetric => (n - 1).toDouble
    var i = 0
    while i < n do
      val x = 2.0 * math.Pi * i.toDouble / denom
      out(i) = a0 - a1 * math.cos(x) + a2 * math.cos(2.0 * x)
      i += 1
    out.result()

  /** SciPy `windows.kaiser`: periodic = length+1 symmetric then drop last. */
  private def kaiser(n: Int, beta: Double, convention: WindowConvention): DVec =
    convention match
      case WindowConvention.Symmetric =>
        kaiserSymmetric(n, beta)
      case WindowConvention.Periodic =>
        val full = kaiserSymmetric(n + 1, beta)
        val out = DVecBuilder.zeros(n)
        var i = 0
        while i < n do
          out(i) = full(i)
          i += 1
        out.result()

  private def kaiserSymmetric(n: Int, beta: Double): DVec =
    val out = DVecBuilder.zeros(n)
    if n == 1 then
      out(0) = 1.0
      out.result()
    else
      val alpha = (n - 1).toDouble / 2.0
      val i0Beta = besselI0(beta)
      var i = 0
      while i < n do
        val t = (i.toDouble - alpha) / alpha
        val r = math.sqrt(math.max(0.0, 1.0 - t * t))
        out(i) = besselI0(beta * r) / i0Beta
        i += 1
      out.result()

  /** Modified Bessel function I₀ via series (SciPy-compatible for window betas). */
  private def besselI0(x: Double): Double =
    val ax = math.abs(x)
    if ax < 3.75 then
      val y = (x / 3.75)
      val y2 = y * y
      1.0 + y2 * (3.5156229 + y2 * (3.0899424 + y2 * (1.2067492 +
        y2 * (0.2659732 + y2 * (0.0360768 + y2 * 0.0045813)))))
    else
      val y = 3.75 / ax
      val expr =
        0.39894228 + y * (0.01328592 + y * (0.00225319 + y * (-0.00157565 +
          y * (0.00916281 + y * (-0.02057706 + y * (0.02635537 +
            y * (-0.01647633 + y * 0.00392377)))))))
      math.exp(ax) / math.sqrt(ax) * expr
