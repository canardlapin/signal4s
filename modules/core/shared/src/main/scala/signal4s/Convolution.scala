package signal4s

import gale.linalg.DVec
import signal4s.internal.{ConvolutionHooks, DirectConvolution}

object Convolution:

  def apply(
      signal: DVec,
      kernel: Kernel,
      region: OutputRegion = OutputRegion.Full,
      method: ConvolutionMethod = ConvolutionMethod.Auto
  ): Either[SignalError, DVec] =
    method match
      case ConvolutionMethod.Direct =>
        Right(direct(signal, kernel, region))
      case ConvolutionMethod.Auto =>
        ConvolutionHooks.runner match
          case Some(run) => run(signal, kernel, region, ConvolutionMethod.Auto)
          case None      => Right(direct(signal, kernel, region))
      case other =>
        ConvolutionHooks.runner match
          case Some(run) => run(signal, kernel, region, other)
          case None =>
            Left(
              SignalError.NumericalFailure(
                "convolution",
                s"method $other requires signal4s-fft (hooks not installed)"
              )
            )

  def apply(
      signal: Signal,
      kernel: Kernel,
      region: OutputRegion,
      method: ConvolutionMethod
  ): Either[SignalError, Signal] =
    apply(signal.samples, kernel, region, method).flatMap { samples =>
      if samples.length == 0 then Left(SignalError.EmptyConvolution)
      else
        val start = outputStart(signal, kernel, region)
        Signal(samples, Sampling(signal.sampleRate, start))
    }

  def apply(
      signal: Signal,
      kernel: Kernel,
      region: OutputRegion
  ): Either[SignalError, Signal] =
    apply(signal, kernel, region, ConvolutionMethod.Auto)

  /** Build a reusable plan (transformed kernel when FFT/OLA is selected). */
  def plan(
      kernel: Kernel,
      inputLength: Int,
      region: OutputRegion = OutputRegion.Full,
      method: ConvolutionMethod = ConvolutionMethod.Auto
  ): Either[SignalError, ConvolutionPlan] =
    ConvolutionHooks.planner match
      case Some(p) => p(kernel, inputLength, region, method)
      case None =>
        method match
          case ConvolutionMethod.Direct | ConvolutionMethod.Auto =>
            Right(DirectConvolutionPlan(kernel, inputLength, region))
          case other =>
            Left(
              SignalError.NumericalFailure(
                "Convolution.plan",
                s"method $other requires signal4s-fft (hooks not installed)"
              )
            )

  /** Circular convolution on \(\mathbb Z/N\mathbb Z\). Not a [[Boundary]] option. */
  def circular(
      signal: DVec,
      kernel: Kernel,
      period: Int
  ): Either[SignalError, DVec] =
    if period <= 0 then Left(SignalError.InvalidPeriod(period))
    else if signal.length != period then
      Left(SignalError.LengthMismatch(period, signal.length))
    else
      ConvolutionHooks.circularRunner match
        case Some(run) => run(signal, kernel, period)
        case None      => Right(DirectConvolution.circular(signal, kernel, period))

  /** Finite batch convolution as a Gale operator with a real adjoint.
    *
    * Supported regions: [[OutputRegion.Full]], [[OutputRegion.Valid]], and
    * [[OutputRegion.Input]] with [[Boundary.Zero]] only. Other boundaries are
    * rejected until their adjoints are implemented.
    */
  def operator(
      kernel: Kernel,
      inputLength: Int,
      region: OutputRegion
  ): Either[SignalError, ConvolutionOperator] =
    if inputLength <= 0 then Left(SignalError.InvalidInputLength(inputLength))
    else
      region match
        case OutputRegion.Full | OutputRegion.Valid | OutputRegion.Input(Boundary.Zero) =>
          val rows = ConvolutionOperator.outputLength(inputLength, kernel, region)
          if rows == 0 then Left(SignalError.EmptyConvolution)
          else Right(ConvolutionOperator.make(kernel, inputLength, region))
        case other =>
          Left(SignalError.UnsupportedOperatorRegion(other))

  private def direct(
      signal: DVec,
      kernel: Kernel,
      region: OutputRegion
  ): DVec =
    region match
      case OutputRegion.Full =>
        DirectConvolution.full(signal, kernel)
      case OutputRegion.Valid =>
        DirectConvolution.valid(signal, kernel)
      case OutputRegion.Input(boundary) =>
        DirectConvolution.inputAligned(signal, kernel, boundary)

  private def outputStart(
      signal: Signal,
      kernel: Kernel,
      region: OutputRegion
  ): Seconds =
    val fs = signal.sampleRate
    region match
      case OutputRegion.Full =>
        DirectConvolution.outputStart(signal.start, kernel, fs)
      case OutputRegion.Valid =>
        val first = DirectConvolution.validFirstSampleIndex(kernel)
        Seconds.unsafe(signal.start.value + first.toDouble / fs.hertz)
      case OutputRegion.Input(_) =>
        signal.start

/** Core-only plan used when FFT hooks are absent. */
private final class DirectConvolutionPlan(
    val kernel: Kernel,
    val inputLength: Int,
    val region: OutputRegion
) extends ConvolutionPlan:
  val selectedMethod: ConvolutionMethod = ConvolutionMethod.Direct

  def apply(signal: DVec): Either[SignalError, DVec] =
    if signal.length != inputLength then
      Left(SignalError.LengthMismatch(inputLength, signal.length))
    else
      Right(
        region match
          case OutputRegion.Full =>
            DirectConvolution.full(signal, kernel)
          case OutputRegion.Valid =>
            DirectConvolution.valid(signal, kernel)
          case OutputRegion.Input(boundary) =>
            DirectConvolution.inputAligned(signal, kernel, boundary)
      )
