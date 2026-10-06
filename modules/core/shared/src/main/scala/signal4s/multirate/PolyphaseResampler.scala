package signal4s.multirate

import gale.linalg.DVec
import signal4s.SignalError

/** Streaming rational resampler: upsample → FIR → downsample.
  *
  * Single-owner, not thread-safe. [[consume]] followed by [[flush]] reproduces
  * batch [[Upfirdn]] on the concatenated input (SciPy `upfirdn` semantics:
  * zero-insertion between samples only — not after the last sample — then FIR
  * transient).
  *
  * `delay` estimates symmetric linear-phase prototype delay in input samples,
  * `(h.length - 1) / (2 * up)`; asymmetric prototypes need their own delay model.
  * Factors are GCD-reduced, as in [[RateRatio]]. Clock/count capacity is Long;
  * exhaustion is refused before changing state.
  */
final class PolyphaseResampler private (
    val ratio: RateRatio,
    val prototype: DVec,
    private val h: DVec,
    private var inputCount: Long,
    private var upIndex: Long,
    private val delayLine: Array[Double],
    private var delayPos: Int,
    private var flushed: Boolean
):
  def up: Int = ratio.up
  def down: Int = ratio.down
  /** Phase of the next high-rate clock modulo the reduced interpolation factor. */
  def phase: Int = (upIndex % up.toLong).toInt
  def delay: Double = (h.length - 1).toDouble / (2.0 * up)
  def samplesConsumed: Long = inputCount
  def isFlushed: Boolean = flushed

  /** Copy delay registers into an owned checkpoint; later calls cannot change it. */
  def snapshot: ResamplerState =
    new ResamplerState(ratio, prototype, DVec.tabulate(delayLine.length)(delayLine(_)),
      delayPos, inputCount, upIndex, flushed)

  /** Restore a compatible stream, including global phase and closed state. */
  def restore(checkpoint: ResamplerState): Either[SignalError, Unit] =
    if checkpoint.ratio != ratio || !samePrototype(checkpoint.prototype) then
      Left(SignalError.NumericalFailure("PolyphaseResampler.restore", "checkpoint ratio/prototype mismatch"))
    else
      var i = 0
      while i < delayLine.length do
        delayLine(i) = checkpoint.registers(i)
        i += 1
      delayPos = checkpoint.delayPosition
      inputCount = checkpoint.samplesConsumed
      upIndex = checkpoint.clock
      flushed = checkpoint.isFlushed
      Right(())

  /** Explicitly start a new segment at clock zero and with empty delay state. */
  def reset(): Unit =
    var i = 0
    while i < delayLine.length do
      delayLine(i) = 0.0
      i += 1
    delayPos = 0
    inputCount = 0L
    upIndex = 0L
    flushed = false

  private def samePrototype(other: DVec): Boolean =
    if other.length != prototype.length then false
    else
      var i = 0
      while i < prototype.length do
        if java.lang.Double.doubleToRawLongBits(other(i)) != java.lang.Double.doubleToRawLongBits(prototype(i)) then
          return false
        i += 1
      true

  def consume(input: DVec): Either[SignalError, DVec] =
    if flushed then
      Left(SignalError.NumericalFailure("PolyphaseResampler", "already flushed"))
    else if !canConsume(input.length) then
      Left(SignalError.NumericalFailure("PolyphaseResampler", "stream clock/count capacity exceeded"))
    else
      val out = scala.collection.mutable.ArrayBuffer.empty[Double]
      var n = 0
      while n < input.length do
        acceptSample(input(n), out)
        n += 1
      Right(DVec.fromSeq(out.toSeq))

  def flush(): Either[SignalError, DVec] =
    if flushed then Right(DVec.zeros(0))
    else if inputCount > 0 && h.length - 1L > Long.MaxValue - upIndex then
      Left(SignalError.NumericalFailure("PolyphaseResampler", "stream clock capacity exceeded during flush"))
    else
      flushed = true
      val out = scala.collection.mutable.ArrayBuffer.empty[Double]
      if inputCount == 0 then Right(DVec.zeros(0))
      else
        // FIR drain only — no up-1 zeros after the final sample (matches SciPy xu length).
        var t = 0
        while t < h.length - 1 do
          clockZero(out)
          t += 1
        Right(DVec.fromSeq(out.toSeq))

  private def canConsume(length: Int): Boolean =
    val clocks =
      if length == 0 then 0L
      else if inputCount == 0 then (length - 1L) * up + 1L
      else length.toLong * up
    length.toLong <= Long.MaxValue - inputCount && clocks <= Long.MaxValue - upIndex

  private def acceptSample(x: Double, out: scala.collection.mutable.ArrayBuffer[Double]): Unit =
    if inputCount > 0 then
      var z = 1
      while z < up do
        clockZero(out)
        z += 1
    inputCount += 1
    pushAndMaybeEmit(x, out)

  private def clockZero(out: scala.collection.mutable.ArrayBuffer[Double]): Unit =
    pushAndMaybeEmit(0.0, out)

  private def pushAndMaybeEmit(
      x: Double,
      out: scala.collection.mutable.ArrayBuffer[Double]
  ): Unit =
    delayLine(delayPos) = x
    var acc = 0.0
    var k = 0
    while k < h.length do
      val idx = Math.floorMod(delayPos - k, delayLine.length)
      acc += h(k) * delayLine(idx)
      k += 1
    delayPos = (delayPos + 1) % delayLine.length
    if upIndex % down.toLong == 0L then out += acc
    upIndex += 1

object PolyphaseResampler:
  def apply(h: DVec, ratio: RateRatio): Either[SignalError, PolyphaseResampler] =
    if h.length == 0 then Left(SignalError.EmptyKernel)
    else
      Right(
        new PolyphaseResampler(
          ratio = ratio,
          prototype = h,
          h = h,
          inputCount = 0L,
          upIndex = 0L,
          delayLine = Array.fill(h.length)(0.0),
          delayPos = 0,
          flushed = false
        )
      )

  def apply(h: DVec, up: Int, down: Int): Either[SignalError, PolyphaseResampler] =
    RateRatio(up, down).flatMap(apply(h, _))
