package signal4s.design

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*
import signal4s.fft.{WindowConvention, WindowSpec}
import signal4s.multirate.{RateRatio, Upfirdn}

/** SciPy-compatible `resample_poly` built on `Upfirdn` + windowed-sinc design. */
object ResamplePoly:

  /** Resample by `up`/`down` (GCD-reduced) with a Kaiser (β=5) lowpass by default.
    *
    * Matches SciPy `resample_poly(..., window=('kaiser', 5.0), padtype='constant')`
    * including filter gain `* up`, pre-pad centering, and output trimming.
    */
  def apply(
      signal: DVec,
      up: Int,
      down: Int,
      window: WindowSpec = WindowSpec.Kaiser(1, 5.0, WindowConvention.Symmetric)
  ): Either[SignalError, DVec] =
    RateRatio(up, down).flatMap { ratio =>
      if signal.length == 0 || ratio.isIdentity then Right(signal.copy)
      else
        val maxRate = math.max(ratio.up, ratio.down)
        val halfLen = 10L * maxRate
        val numTaps = 2L * halfLen + 1L
        val nPrePad = ratio.down.toLong - halfLen % ratio.down
        val nPreRemove = (halfLen + nPrePad) / ratio.down
        val product = signal.length.toLong * ratio.up
        val nOut = product / ratio.down + (if product % ratio.down == 0 then 0L else 1L)
        val needed = nOut + nPreRemove
        // Solve the full upfirdn length inequality directly, avoiding a padding loop.
        val minimumFilterLength = (needed - 1L) * ratio.down - (signal.length - 1L) * ratio.up + 1L
        val nPostPad = math.max(0L, minimumFilterLength - numTaps - nPrePad)
        val hLen = numTaps + nPrePad + nPostPad
        if numTaps > Int.MaxValue || needed > Int.MaxValue || hLen > Int.MaxValue then
          Left(SignalError.NumericalFailure("ResamplePoly", "prototype or output length exceeds Int capacity"))
        else
          // firwin cutoff is relative to Nyquist; at fs=2, Nyquist=1.
          val fs = SampleRate.unsafe(2.0)
          val cutoff = Frequency.unsafe(1.0 / maxRate)
          val win = withLength(window, numTaps.toInt)
          FirDesign.lowPass(numTaps.toInt, cutoff, fs, win).flatMap { designed =>
            val h0 = DVec.tabulate(designed.fir.taps.length)(i => designed.fir.taps(i) * ratio.up)
            val h = padFilter(h0, nPrePad.toInt, nPostPad.toInt)
            Upfirdn(h, signal, ratio).map { y =>
              val until = math.min(y.length.toLong, needed).toInt
              if nPreRemove >= until then DVec.zeros(0)
              else y.slice(nPreRemove.toInt, until).copy
            }
          }
    }

  private def padFilter(h: DVec, pre: Int, post: Int): DVec =
    val out = DVecBuilder.zeros(h.length + pre + post)
    var i = 0
    while i < h.length do
      out(pre + i) = h(i)
      i += 1
    out.result()

  private def withLength(spec: WindowSpec, n: Int): WindowSpec =
    spec match
      case WindowSpec.Rectangular(_, c)  => WindowSpec.Rectangular(n, c)
      case WindowSpec.Hann(_, c)         => WindowSpec.Hann(n, c)
      case WindowSpec.Hamming(_, c)      => WindowSpec.Hamming(n, c)
      case WindowSpec.Blackman(_, c)     => WindowSpec.Blackman(n, c)
      case WindowSpec.Kaiser(_, beta, c) => WindowSpec.Kaiser(n, beta, c)
