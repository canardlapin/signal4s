package signal4s.filter

import gale.linalg.{DVec, MutableDVec}
import signal4s.*

/** Causal FIR system described by a [[Kernel]] with `zeroLagIndex == 0`.
  *
  * Streaming runners are single-owner and not thread-safe. Algorithmic (group)
  * delay for a length-`M` symmetric linear-phase FIR is `(M-1)/2` samples;
  * this type does not assume symmetry or assign a universal group delay.
  * Causal support extends from lag 0 to `M-1`. State length is
  * `M-1` (DF-II transposed delays).
  */
final case class Fir private (kernel: Kernel):
  def taps: DVec = kernel.taps
  def length: Int = kernel.length
  def stateLength: Int = math.max(0, kernel.length - 1)

  def newRunner(): FirRunner =
    FirRunner(this)

  def newRunner(initial: FilterState): Either[SignalError, FirRunner] =
    FirRunner(this, initial)

  /** Batch convenience: filter from rest. */
  def process(samples: DVec): Either[SignalError, DVec] =
    newRunner().process(samples)

  def process(signal: Signal): Either[SignalError, Signal] =
    process(signal.samples).flatMap(out => Signal(out, signal.sampling))

object Fir:
  def apply(kernel: Kernel): Either[SignalError, Fir] =
    if kernel.zeroLagIndex != 0 then
      Left(
        SignalError.NumericalFailure(
          "Fir",
          s"FIR runners require a causal kernel (zeroLagIndex=0), got ${kernel.zeroLagIndex}"
        )
      )
    else Right(new Fir(kernel))

  def causal(taps: DVec): Either[SignalError, Fir] =
    Kernel.causal(taps).flatMap(apply)

/** Single-owner FIR runner (DF-II transposed). Not thread-safe. */
final class FirRunner private (
    val fir: Fir,
    private val state: Array[Double]
):
  def stateLength: Int = fir.stateLength

  def process(input: DVec): Either[SignalError, DVec] =
    val out = MutableDVec.zeros(input.length)
    processInto(input, out).map(_ => out.toVec)

  def processInto(input: DVec, output: MutableDVec): Either[SignalError, Unit] =
    if output.length != input.length then
      Left(SignalError.LengthMismatch(input.length, output.length))
    else
      val b = fir.taps
      Df2Transposed.filterInto(
        input,
        output,
        x => Df2Transposed.transferStep(b, FirRunner.AOne, x, state)
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

object FirRunner:
  private val AOne: DVec = DVec.tabulate(1)(_ => 1.0)

  private[filter] def apply(fir: Fir): FirRunner =
    new FirRunner(fir, Array.fill(fir.stateLength)(0.0))

  private[filter] def apply(
      fir: Fir,
      initial: FilterState
  ): Either[SignalError, FirRunner] =
    if initial.length != fir.stateLength then
      Left(SignalError.InvalidFilterState(fir.stateLength, initial.length))
    else
      val st = Array.tabulate(fir.stateLength)(i => initial.values(i))
      Right(new FirRunner(fir, st))
