package signal4s.design

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*
import signal4s.fft.{WindowConvention, WindowSpec}
import signal4s.multirate.{RateRatio, Upfirdn}

/** SciPy-compatible `resample_poly` built on [[Upfirdn]] + windowed-sinc design. */
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
      if ratio.isIdentity then Right(signal.copy)
      else
        val maxRate = math.max(ratio.up, ratio.down)
        val halfLen = 10 * maxRate
        val numTaps = 2 * halfLen + 1
        // firwin cutoff is relative to Nyquist; use Frequency at fs=2 so Nyquist=1
        val fs = SampleRate.unsafe(2.0)
        val cutoff = Frequency.unsafe(1.0 / maxRate) // cycles/sample * Nyquist=1 → Hz at fs=2
        val win = withLength(window, numTaps)
        FirDesign.lowPass(numTaps, cutoff, fs, win).flatMap { designed =>
          // Scale by up (SciPy)
          val h0 = DVec.tabulate(designed.fir.taps.length)(i => designed.fir.taps(i) * ratio.up)
          val nPrePad = ratio.down - halfLen % ratio.down
          val nOut = outputLen(signal.length, ratio.up, ratio.down)
          val nPreRemove = (halfLen + nPrePad) / ratio.down
          var nPostPad = 0
          var hLen = h0.length + nPrePad + nPostPad
          while Upfirdn.outputLength(hLen, signal.length, ratio.up, ratio.down) < nOut + nPreRemove do
            nPostPad += 1
            hLen = h0.length + nPrePad + nPostPad
          val h = padFilter(h0, nPrePad, nPostPad)
          Upfirdn(h, signal, ratio).map { y =>
            val until = math.min(y.length, nPreRemove + nOut)
            if nPreRemove >= until then DVec.zeros(0)
            else y.slice(nPreRemove, until).copy
          }
        }
    }

  private def outputLen(nIn: Int, up: Int, down: Int): Int =
    val nOut = nIn.toLong * up
    (nOut / down + (if nOut % down != 0 then 1 else 0)).toInt

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
