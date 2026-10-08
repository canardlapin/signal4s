package signal4s.design

import gale.linalg.DVec
import signal4s.*
import signal4s.fft.{WindowConvention, WindowSpec}
import signal4s.multirate.RateRatio

/** SciPy-compatible centered `resample_poly`, sharing the bounded plan's numeric core. */
object ResamplePoly:
  def apply(signal: DVec,up: Int,down: Int,
      window: WindowSpec=WindowSpec.Kaiser(1,5.0,WindowConvention.Symmetric)): Either[SignalError,DVec] =
    RateRatio(up,down).flatMap: ratio =>
      // Preserve empty-input admission without constructing a needless prototype.
      if signal.length==0 then Right(signal.copy)
      else ResamplePolyPlan(ratio.up,ratio.down,window).flatMap: plan =>
        plan.outputLength(signal.length).flatMap: count =>
          if count>Int.MaxValue then Left(SignalError.NumericalFailure("ResamplePoly","output length exceeds Int capacity"))
          else plan.processWindow(signal,0L,signal.length,0L,count.toInt)
