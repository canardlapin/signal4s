package signal4s.filter

import gale.linalg.{DVec, MutableDVec}
import signal4s.*

/** Digital transfer function in negative powers of \(z\):
  * \( A(z^{-1}) Y = B(z^{-1}) X \) with `a(0)` normalized to 1.
  */
final case class DigitalTransferFunction private (
    feedForward: DVec,
    feedback: DVec
):
  def b: DVec = feedForward
  def a: DVec = feedback
  def stateLength: Int = math.max(feedForward.length, feedback.length) - 1

  def newRunner(): TransferFunctionRunner =
    TransferFunctionRunner(this)

  def newRunner(initial: FilterState): Either[SignalError, TransferFunctionRunner] =
    TransferFunctionRunner(this, initial)

  def process(samples: DVec): Either[SignalError, DVec] =
    newRunner().process(samples)

  def process(signal: Signal): Either[SignalError, Signal] =
    process(signal.samples).flatMap(out => Signal(out, signal.sampling))

object DigitalTransferFunction:
  def apply(
      feedForward: DVec,
      feedback: DVec
  ): Either[SignalError, DigitalTransferFunction] =
    if feedForward.length == 0 then
      Left(SignalError.InvalidFilterCoefficient("b", Double.NaN))
    else if feedback.length == 0 then
      Left(SignalError.InvalidFilterCoefficient("a", Double.NaN))
    else if !feedback(0).isFinite || feedback(0) == 0.0 then
      Left(SignalError.InvalidFilterCoefficient("a0", feedback(0)))
    else
      val a0 = feedback(0)
      val b =
        if a0 == 1.0 then feedForward
        else DVec.tabulate(feedForward.length)(i => feedForward(i) / a0)
      val a =
        if a0 == 1.0 then feedback
        else DVec.tabulate(feedback.length)(i => feedback(i) / a0)
      var i = 0
      while i < b.length do
        if !b(i).isFinite then
          return Left(SignalError.InvalidFilterCoefficient(s"b$i", b(i)))
        i += 1
      i = 0
      while i < a.length do
        if !a(i).isFinite then
          return Left(SignalError.InvalidFilterCoefficient(s"a$i", a(i)))
        i += 1
      Right(new DigitalTransferFunction(b, a))

final class TransferFunctionRunner private (
    val system: DigitalTransferFunction,
    private val state: Array[Double]
):
  def stateLength: Int = system.stateLength

  def process(input: DVec): Either[SignalError, DVec] =
    val out = MutableDVec.zeros(input.length)
    processInto(input, out).map(_ => out.toVec)

  def processInto(input: DVec, output: MutableDVec): Either[SignalError, Unit] =
    if output.length != input.length then
      Left(SignalError.LengthMismatch(input.length, output.length))
    else
      Df2Transposed.filterInto(
        input,
        output,
        x => Df2Transposed.transferStep(system.b, system.a, x, state)
      )
      Right(())

  def reset(): Unit =
    var i = 0
    while i < state.length do
      state(i) = 0.0
      i += 1

  def snapshot: FilterState =
    FilterState.unsafe(DVec.tabulate(state.length)(state(_)))

  def restore(snapshot: FilterState): Either[SignalError, Unit] =
    if snapshot.length != state.length then
      Left(SignalError.InvalidFilterState(state.length, snapshot.length))
    else
      var i = 0
      while i < state.length do
        state(i) = snapshot.values(i)
        i += 1
      Right(())

object TransferFunctionRunner:
  private[filter] def apply(system: DigitalTransferFunction): TransferFunctionRunner =
    new TransferFunctionRunner(system, Array.fill(system.stateLength)(0.0))

  private[filter] def apply(
      system: DigitalTransferFunction,
      initial: FilterState
  ): Either[SignalError, TransferFunctionRunner] =
    if initial.length != system.stateLength then
      Left(SignalError.InvalidFilterState(system.stateLength, initial.length))
    else
      Right(
        new TransferFunctionRunner(
          system,
          Array.tabulate(system.stateLength)(i => initial.values(i))
        )
      )
