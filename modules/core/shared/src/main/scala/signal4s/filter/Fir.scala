package signal4s.filter

import gale.linalg.{DVec, MutableDVec}
import signal4s.*

/** Primitive numerical payloads, excluding caller input/destination and object/GC/RSS overhead. */
final case class FirResources(coefficientBytes: BigInt,stateBytes: BigInt,snapshotBytes: BigInt):
  def processIntoScratchBytes(samples: Int): Either[SignalError,BigInt] =
    if samples<0 then Left(SignalError.InvalidInputLength(samples)) else Right(stateBytes+BigInt(samples)*8)
  def ownedProcessAdditionalBytes(samples: Int): Either[SignalError,BigInt] =
    if samples<0 then Left(SignalError.InvalidInputLength(samples)) else Right(stateBytes+BigInt(samples)*24)

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
  def resources: FirResources=FirResources((BigInt(length)+1)*8,BigInt(stateLength)*8,BigInt(stateLength)*8)

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
  private[filter] def finite(values: DVec): Boolean =
    var i=0
    while i<values.length do
      if !values(i).isFinite then return false
      i+=1
    true
  def apply(kernel: Kernel): Either[SignalError, Fir] =
    if kernel.length<=0 then Left(SignalError.EmptyKernel)
    else if !Fir.finite(kernel.taps) then Left(SignalError.NumericalFailure("Fir","finite coefficients required"))
    else if kernel.zeroLagIndex != 0 then
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
      if b.length<=0 || !Fir.finite(b) || !state.forall(_.isFinite) || !Fir.finite(input) then
        Left(SignalError.NumericalFailure("FirRunner","finite coefficients, state and input required"))
      else
        val next=state.clone()
        val staged=Array.ofDim[Double](input.length)
        var i=0
        var failed=false
        while i<input.length && !failed do
          val y=Df2Transposed.transferStep(b,FirRunner.AOne,input(i),next)
          if !y.isFinite || !next.forall(_.isFinite) then failed=true else staged(i)=y
          i+=1
        if failed then Left(SignalError.NumericalFailure("FirRunner","output/delay arithmetic exceeds finite capacity"))
        else
          i=0
          while i<input.length do {output(i)=staged(i);i+=1}
          i=0
          while i<state.length do {state(i)=next(i);i+=1}
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
    else if !Fir.finite(snapshot.values) then
      Left(SignalError.NumericalFailure("FirRunner.restore","finite delay state required"))
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
    else if !Fir.finite(initial.values) then Left(SignalError.NumericalFailure("FirRunner","finite initial state required"))
    else
      val st = Array.tabulate(fir.stateLength)(i => initial.values(i))
      Right(new FirRunner(fir, st))
