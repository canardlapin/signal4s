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
  * @param delay linear-phase group-delay estimate in input samples,
  *              \((\texttt{h.length}-1)/2\).
  */
final class PolyphaseResampler private (
    val ratio: RateRatio,
    val prototype: DVec,
    private val h: DVec,
    private var inputCount: Int,
    private var upIndex: Int,
    private val delayLine: Array[Double],
    private var delayPos: Int,
    private var flushed: Boolean
):
  def up: Int = ratio.up
  def down: Int = ratio.down
  def phase: Int = Math.floorMod(upIndex, up)
  def delay: Double = (h.length - 1).toDouble / 2.0
  def samplesConsumed: Int = inputCount
  def isFlushed: Boolean = flushed

  def consume(input: DVec): Either[SignalError, DVec] =
    if flushed then
      Left(SignalError.NumericalFailure("PolyphaseResampler", "already flushed"))
    else
      val out = scala.collection.mutable.ArrayBuffer.empty[Double]
      var n = 0
      while n < input.length do
        acceptSample(input(n), out)
        n += 1
      Right(DVec.fromSeq(out.toSeq))

  def flush(): Either[SignalError, DVec] =
    if flushed then Right(DVec.zeros(0))
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
    if Math.floorMod(upIndex, down) == 0 then out += acc
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
          inputCount = 0,
          upIndex = 0,
          delayLine = Array.fill(h.length)(0.0),
          delayPos = 0,
          flushed = false
        )
      )

  def apply(h: DVec, up: Int, down: Int): Either[SignalError, PolyphaseResampler] =
    RateRatio(up, down).flatMap(apply(h, _))
