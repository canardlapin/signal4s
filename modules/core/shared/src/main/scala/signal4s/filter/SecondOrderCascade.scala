package signal4s.filter

import gale.linalg.{DVec, MutableDVec}
import signal4s.*

/** Cascade of normalized biquads with an overall gain.
  *
  * Execution uses direct-form II transposed sections. Runners are single-owner
  * and not thread-safe. State length is `2 * sections.length`. There is no
  * end-of-stream flush beyond consuming the provided samples.
  */
final case class SecondOrderCascade private (
    gain: Double,
    sections: IArray[Biquad]
):
  def sectionCount: Int = sections.length
  def stateLength: Int = 2 * sections.length

  def newRunner(): SosRunner =
    SosRunner(this)

  def newRunner(initial: FilterState): Either[SignalError, SosRunner] =
    SosRunner(this, initial)

  def process(samples: DVec): Either[SignalError, DVec] =
    newRunner().process(samples)

  def process(signal: Signal): Either[SignalError, Signal] =
    process(signal.samples).flatMap(out => Signal(out, signal.sampling))

object SecondOrderCascade:
  def apply(
      gain: Double,
      sections: IArray[Biquad]
  ): Either[SignalError, SecondOrderCascade] =
    if !gain.isFinite then Left(SignalError.InvalidFilterCoefficient("gain", gain))
    else Right(new SecondOrderCascade(gain, sections))

  def of(sections: Biquad*): Either[SignalError, SecondOrderCascade] =
    apply(1.0, IArray.from(sections))

  /** Build from SciPy SOS matrix rows of length 6; overall gain stays 1. */
  def fromSosMatrix(rows: IArray[IArray[Double]]): Either[SignalError, SecondOrderCascade] =
    val built = IArray.newBuilder[Biquad]
    var i = 0
    while i < rows.length do
      Biquad.fromSosRow(rows(i)) match
        case Left(err) => return Left(err)
        case Right(bq) => built += bq
      i += 1
    apply(1.0, built.result())

final class SosRunner private (
    val cascade: SecondOrderCascade,
    private val state: Array[Double]
):
  def stateLength: Int = cascade.stateLength

  def process(input: DVec): Either[SignalError, DVec] =
    val out = MutableDVec.zeros(input.length)
    processInto(input, out).map(_ => out.toVec)

  def processInto(input: DVec, output: MutableDVec): Either[SignalError, Unit] =
    if output.length != input.length then
      Left(SignalError.LengthMismatch(input.length, output.length))
    else
      var n = 0
      while n < input.length do
        var v = cascade.gain * input(n)
        var s = 0
        while s < cascade.sections.length do
          val base = 2 * s
          // local 2-state view
          val d0 = state(base)
          val d1 = state(base + 1)
          val sec = cascade.sections(s)
          val y = sec.b0 * v + d0
          state(base) = sec.b1 * v - sec.a1 * y + d1
          state(base + 1) = sec.b2 * v - sec.a2 * y
          v = y
          s += 1
        output(n) = v
        n += 1
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

object SosRunner:
  private[filter] def apply(cascade: SecondOrderCascade): SosRunner =
    new SosRunner(cascade, Array.fill(cascade.stateLength)(0.0))

  private[filter] def apply(
      cascade: SecondOrderCascade,
      initial: FilterState
  ): Either[SignalError, SosRunner] =
    if initial.length != cascade.stateLength then
      Left(SignalError.InvalidFilterState(cascade.stateLength, initial.length))
    else
      Right(
        new SosRunner(
          cascade,
          Array.tabulate(cascade.stateLength)(i => initial.values(i))
        )
      )
