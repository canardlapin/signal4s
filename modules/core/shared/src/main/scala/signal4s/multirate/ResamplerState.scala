package signal4s.multirate

import gale.linalg.DVec

/** Owned immutable checkpoint of a single resampler stream.
  *
  * Registers are an independent snapshot. A checkpoint is specific to the exact
  * prototype and reduced ratio, and includes global clocks and flush state.
  * It carries no acquisition-segment identity or serialized-format promise.
  */
final class ResamplerState private[multirate] (
    val ratio: RateRatio,
    private[multirate] val prototype: DVec,
    val registers: DVec,
    private[multirate] val delayPosition: Int,
    val samplesConsumed: Long,
    private[multirate] val clock: Long,
    val isFlushed: Boolean
)
