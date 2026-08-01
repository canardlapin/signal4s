package signal4s.internal

import gale.linalg.DVec
import signal4s.*

/** Optional FFT/OLA backend installed by `signal4s-fft` without a core→fft dependency. */
private[signal4s] object ConvolutionHooks:

  type Runner =
    (DVec, Kernel, OutputRegion, ConvolutionMethod) => Either[SignalError, DVec]

  type Planner =
    (Kernel, Int, OutputRegion, ConvolutionMethod) => Either[SignalError, ConvolutionPlan]
  type CircularRunner = (DVec, Kernel, Int) => Either[SignalError, DVec]

  @volatile var runner: Option[Runner] = None
  @volatile var planner: Option[Planner] = None
  @volatile var circularRunner: Option[CircularRunner] = None

  def install(run: Runner, plan: Planner, circular: CircularRunner): Unit =
    runner = Some(run)
    planner = Some(plan)
    circularRunner = Some(circular)
