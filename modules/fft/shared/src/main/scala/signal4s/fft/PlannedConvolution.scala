package signal4s.fft

import gale.linalg.DVec
import signal4s.*

/** [[ConvolutionPlan]] that reuses a transformed kernel spectrum when using FFT/OLA.
  *
  * FFT and overlap-add plans retain mutable FFT and output-accumulation buffers.
  * [[apply]] is therefore single-owner / not thread-safe for those methods (same
  * discipline as filter runners). Direct plans remain pure.
  */
final class PlannedConvolution private (
    val kernel: Kernel,
    val inputLength: Int,
    val region: OutputRegion,
    val selectedMethod: ConvolutionMethod,
    private val kernelSpectrumArrays: Option[(Array[Double], Array[Double])],
    private val workRe: Array[Double],
    private val workIm: Array[Double],
    private val packRe: Array[Double],
    private val packIm: Array[Double],
    private val outBuf: Array[Double],
    private val accumulationBuf: Array[Double]
) extends ConvolutionPlan:

  def apply(signal: DVec): Either[SignalError, DVec] =
    if signal.length != inputLength then
      Left(SignalError.LengthMismatch(inputLength, signal.length))
    else
      selectedMethod match
        case ConvolutionMethod.Direct =>
          Right(PlannedConvolution.direct(signal, kernel, region))
        case ConvolutionMethod.Fft =>
          kernelSpectrumArrays match
            case None =>
              Left(SignalError.NumericalFailure("ConvolutionPlan", "missing kernel spectrum"))
            case Some((kr, ki)) =>
              val outLen = signal.length + kernel.length - 1
              FftConvolution.regionBounds(signal.length, kernel, region).flatMap {
                case (_, 0) => Right(DVec.zeros(0))
                case (offset, length) =>
                  FftConvolution
                    .fullWithKernelSpectrumArraysIntoRegion(
                      signal,
                      kr,
                      ki,
                      outLen,
                      offset,
                      length,
                      workRe,
                      workIm,
                      packRe,
                      packIm,
                      outBuf
                    )
                    .map(_ => FftConvolution.adoptFull(outBuf, length))
              }
        case ConvolutionMethod.OverlapAdd(block) =>
          kernelSpectrumArrays match
            case None =>
              Left(SignalError.NumericalFailure("ConvolutionPlan", "missing kernel spectrum"))
            case Some((kr, ki)) =>
              FftConvolution.regionBounds(signal.length, kernel, region).flatMap {
                case (_, 0) => Right(DVec.zeros(0))
                case (offset, length) =>
                  OverlapAddConvolution
                    .fullWithKernelSpectrumArraysIntoRegion(
                      signal,
                      kr,
                      ki,
                      kernel.length,
                      block,
                      signal.length + kernel.length - 1,
                      offset,
                      length,
                      workRe,
                      workIm,
                      packRe,
                      packIm,
                      outBuf,
                      accumulationBuf
                    )
                    .map(_ => FftConvolution.adoptFull(accumulationBuf, length))
              }
        case ConvolutionMethod.Auto =>
          Left(SignalError.NumericalFailure("ConvolutionPlan", "Auto must be resolved at plan time"))

object PlannedConvolution:

  def make(
      kernel: Kernel,
      inputLength: Int,
      region: OutputRegion,
      method: ConvolutionMethod
  ): Either[SignalError, PlannedConvolution] =
    if inputLength <= 0 then Left(SignalError.InvalidInputLength(inputLength))
    else if !supports(region, method) then
      method match
        case ConvolutionMethod.Auto | ConvolutionMethod.Direct =>
          Right(directPlan(kernel, inputLength, region))
        case other =>
          Left(
            SignalError.NumericalFailure(
              "ConvolutionPlan",
              s"$other unsupported for region $region"
            )
          )
    else
      val resolved = resolve(method, inputLength, kernel.length, region)
      resolved match
        case ConvolutionMethod.Direct =>
          Right(directPlan(kernel, inputLength, region))
        case ConvolutionMethod.Fft =>
          val outLen = inputLength + kernel.length - 1
          val nfft = FftConvolution.fftLength(outLen)
          for
            bounds <- FftConvolution.regionBounds(inputLength, kernel, region)
            half <- FftConvolution.transformKernelHalfArrays(kernel, nfft)
          yield
            val packN =
              if signal4s.fft.internal.FastFftLength.isPowerOfTwo(nfft) then nfft / 2 else 0
            new PlannedConvolution(
              kernel,
              inputLength,
              region,
              ConvolutionMethod.Fft,
              Some(half),
              new Array[Double](nfft),
              new Array[Double](nfft),
              new Array[Double](packN),
              new Array[Double](packN),
              new Array[Double](bounds._2),
              Array.empty
            )
        case ConvolutionMethod.OverlapAdd(block) =>
          val nfft = FftConvolution.fftLength(block + kernel.length - 1)
          for
            bounds <- FftConvolution.regionBounds(inputLength, kernel, region)
            half <- FftConvolution.transformKernelHalfArrays(kernel, nfft)
          yield
            val packN =
              if signal4s.fft.internal.FastFftLength.isPowerOfTwo(nfft) then nfft / 2 else 0
            new PlannedConvolution(
              kernel,
              inputLength,
              region,
              ConvolutionMethod.OverlapAdd(block),
              Some(half),
              new Array[Double](nfft),
              new Array[Double](nfft),
              new Array[Double](packN),
              new Array[Double](packN),
              new Array[Double](nfft),
              new Array[Double](bounds._2)
            )
        case ConvolutionMethod.Auto =>
          Left(SignalError.NumericalFailure("ConvolutionPlan", "internal: unresolved Auto"))

  private def directPlan(
      kernel: Kernel,
      inputLength: Int,
      region: OutputRegion
  ): PlannedConvolution =
    new PlannedConvolution(
      kernel,
      inputLength,
      region,
      ConvolutionMethod.Direct,
      None,
      Array.empty,
      Array.empty,
      Array.empty,
      Array.empty,
      Array.empty,
      Array.empty
    )

  private def resolve(
      method: ConvolutionMethod,
      n: Int,
      m: Int,
      region: OutputRegion
  ): ConvolutionMethod =
    method match
      case ConvolutionMethod.Auto =>
        if fftCapable(region) then AutoCostModel.select(n, m)
        else ConvolutionMethod.Direct
      case ConvolutionMethod.OverlapAdd(0) =>
        ConvolutionMethod.OverlapAdd(AutoCostModel.defaultBlock(n, m))
      case other => other

  private def supports(region: OutputRegion, method: ConvolutionMethod): Boolean =
    method match
      case ConvolutionMethod.Direct => true
      case ConvolutionMethod.Auto   => true
      case _                        => fftCapable(region)

  private def fftCapable(region: OutputRegion): Boolean =
    region match
      case OutputRegion.Full | OutputRegion.Valid | OutputRegion.Input(Boundary.Zero) =>
        true
      case _ => false

  private[fft] def direct(signal: DVec, kernel: Kernel, region: OutputRegion): DVec =
    import signal4s.internal.DirectConvolution
    region match
      case OutputRegion.Full =>
        DirectConvolution.full(signal, kernel)
      case OutputRegion.Valid =>
        DirectConvolution.valid(signal, kernel)
      case OutputRegion.Input(boundary) =>
        DirectConvolution.inputAligned(signal, kernel, boundary)
