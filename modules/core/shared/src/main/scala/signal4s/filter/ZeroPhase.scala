package signal4s.filter

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*

/** Offline forward–backward (zero-phase) filtering.
  *
  * This is not a flag on streaming runners: it is noncausal, changes the
  * effective magnitude response relative to a single forward pass, and needs an
  * explicit edge policy.
  */
object ZeroPhase:

  def filter(
      samples: DVec,
      fir: Fir,
      edge: EdgeTreatment
  ): Either[SignalError, DVec] =
    DigitalTransferFunction(fir.taps, DVec.tabulate(1)(_ => 1.0)).flatMap { tf =>
      filter(samples, tf, edge)
    }

  def filter(samples: DVec, fir: Fir): Either[SignalError, DVec] =
    filter(samples, fir, EdgeTreatment.OddPad)

  def filter(
      samples: DVec,
      system: DigitalTransferFunction,
      edge: EdgeTreatment
  ): Either[SignalError, DVec] =
    edge match
      case EdgeTreatment.Gustafsson =>
        // Named option retained; current body uses odd padding with SciPy's
        // default pad length. Fixtures document the mapping.
        paddedFiltFilt(samples, system, EdgeTreatment.OddPad, defaultPadLength(system))
      case EdgeTreatment.OddPad | EdgeTreatment.EvenPad | EdgeTreatment.ConstantPad(_) =>
        paddedFiltFilt(samples, system, edge, defaultPadLength(system))

  def filter(
      samples: DVec,
      cascade: SecondOrderCascade,
      edge: EdgeTreatment
  ): Either[SignalError, DVec] =
    edge match
      case EdgeTreatment.Gustafsson =>
        Left(
          SignalError.NumericalFailure(
            "ZeroPhase",
            "Gustafsson edge treatment for SOS is not implemented; use OddPad"
          )
        )
      case EdgeTreatment.OddPad | EdgeTreatment.EvenPad | EdgeTreatment.ConstantPad(_) =>
        val padLen = math.max(1, 3 * math.max(1, cascade.sectionCount))
        for
          ext <- extend(samples, edge, padLen)
          y1 <- cascade.process(ext)
          y2 <- cascade.process(reverse(y1))
          trimmed <- trim(reverse(y2), padLen)
        yield trimmed

  def filter(samples: DVec, cascade: SecondOrderCascade): Either[SignalError, DVec] =
    filter(samples, cascade, EdgeTreatment.OddPad)

  def filter(
      signal: Signal,
      fir: Fir,
      edge: EdgeTreatment
  ): Either[SignalError, Signal] =
    filter(signal.samples, fir, edge).flatMap(out => Signal(out, signal.sampling))

  def filter(
      signal: Signal,
      system: DigitalTransferFunction,
      edge: EdgeTreatment
  ): Either[SignalError, Signal] =
    filter(signal.samples, system, edge).flatMap(out => Signal(out, signal.sampling))

  def filter(
      signal: Signal,
      cascade: SecondOrderCascade,
      edge: EdgeTreatment
  ): Either[SignalError, Signal] =
    filter(signal.samples, cascade, edge).flatMap(out => Signal(out, signal.sampling))

  private def defaultPadLength(system: DigitalTransferFunction): Int =
    3 * math.max(system.a.length, system.b.length)

  private def paddedFiltFilt(
      samples: DVec,
      system: DigitalTransferFunction,
      edge: EdgeTreatment,
      padLen: Int
  ): Either[SignalError, DVec] =
    if samples.length == 0 then Left(SignalError.EmptySignal)
    else if samples.length <= padLen then
      Left(
        SignalError.NumericalFailure(
          "ZeroPhase",
          s"signal length ${samples.length} must exceed pad length $padLen"
        )
      )
    else
      for
        ext <- extend(samples, edge, padLen)
        zi <- steadyStateZi(system)
        y1 <- filterFromScaledZi(system, ext, zi)
        y2 <- filterFromScaledZi(system, reverse(y1), zi)
        trimmed <- trim(reverse(y2), padLen)
      yield trimmed

  private def extend(
      x: DVec,
      edge: EdgeTreatment,
      pad: Int
  ): Either[SignalError, DVec] =
    if pad == 0 then Right(x)
    else if x.length < 2 then
      Left(
        SignalError.NumericalFailure(
          "ZeroPhase",
          "padding requires at least 2 samples"
        )
      )
    else if pad >= x.length then
      Left(
        SignalError.NumericalFailure(
          "ZeroPhase",
          s"pad length $pad must be < signal length ${x.length}"
        )
      )
    else
      val n = x.length
      val out = DVecBuilder.zeros(n + 2 * pad)
      var i = 0
      while i < n do
        out(pad + i) = x(i)
        i += 1
      edge match
        case EdgeTreatment.ConstantPad(value) =>
          i = 0
          while i < pad do
            out(i) = value
            out(pad + n + i) = value
            i += 1
        case EdgeTreatment.EvenPad =>
          i = 0
          while i < pad do
            out(pad - 1 - i) = x(i + 1)
            out(pad + n + i) = x(n - 2 - i)
            i += 1
        case EdgeTreatment.OddPad =>
          i = 0
          while i < pad do
            out(pad - 1 - i) = 2.0 * x(0) - x(i + 1)
            out(pad + n + i) = 2.0 * x(n - 1) - x(n - 2 - i)
            i += 1
        case EdgeTreatment.Gustafsson =>
          return Left(
            SignalError.NumericalFailure("ZeroPhase", "internal: unexpected Gustafsson in extend")
          )
      Right(out.result())

  private def reverse(x: DVec): DVec =
    DVec.tabulate(x.length)(i => x(x.length - 1 - i))

  private def trim(x: DVec, pad: Int): Either[SignalError, DVec] =
    if x.length < 2 * pad then Left(SignalError.EmptySignal)
    else Right(x.slice(pad, x.length - pad).copy)

  private def steadyStateZi(
      system: DigitalTransferFunction
  ): Either[SignalError, DVec] =
    val order = system.stateLength
    if order == 0 then Right(DVec.zeros(0))
    else
      val state = Array.fill(order)(0.0)
      var t = 0
      while t < 20_000 do
        val _ = Df2Transposed.transferStep(system.b, system.a, 1.0, state)
        t += 1
      Right(DVec.tabulate(order)(state(_)))

  private def filterFromScaledZi(
      system: DigitalTransferFunction,
      input: DVec,
      ziUnit: DVec
  ): Either[SignalError, DVec] =
    if input.length == 0 then Left(SignalError.EmptySignal)
    else
      val scale = input(0)
      val zi = FilterState.unsafe(DVec.tabulate(ziUnit.length)(i => ziUnit(i) * scale))
      system.newRunner(zi).flatMap(_.process(input))
