package signal4s.multirate

import gale.linalg.DVec
import signal4s.*

/** Convenience rate-change ops that compile to [[Upfirdn]] (no second numeric core). */
object Resample:

  /** Rational resample with an explicit FIR prototype (SciPy `upfirdn` path). */
  def apply(signal: DVec, h: DVec, ratio: RateRatio): Either[SignalError, DVec] =
    Upfirdn(h, signal, ratio)

  def apply(signal: DVec, h: DVec, up: Int, down: Int): Either[SignalError, DVec] =
    Upfirdn(h, signal, up, down)

  /** Pure integer decimation with FIR (downsample after filter). */
  def decimate(signal: DVec, h: DVec, factor: Int): Either[SignalError, DVec] =
    RateRatio(1, factor).flatMap(r => Upfirdn(h, signal, r))

  /** Pure integer interpolation with FIR (upsample before filter). */
  def interpolate(signal: DVec, h: DVec, factor: Int): Either[SignalError, DVec] =
    RateRatio(factor, 1).flatMap(r => Upfirdn(h, signal, r))

  /** Identity rate change — still applies FIR (documents delay/filter effect). */
  def identityFiltered(signal: DVec, h: DVec): Either[SignalError, DVec] =
    Upfirdn(h, signal, RateRatio.identity)
