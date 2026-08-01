package signal4s.fft

import gale.linalg.DVec
import signal4s.*

/** Overlap-add linear convolution with FFT block processing. */
private[fft] object OverlapAddConvolution:

  def full(
      signal: DVec,
      kernel: Kernel,
      blockLength: Int
  ): Either[SignalError, DVec] =
    if blockLength <= 0 then Left(SignalError.InvalidInputLength(blockLength))
    else
      val m = kernel.length
      val nfft = FftConvolution.fftLength(blockLength + m - 1)
      FftConvolution.transformKernelHalfArrays(kernel, nfft).flatMap { case (kr, ki) =>
        val packN =
          if signal4s.fft.internal.FastFftLength.isPowerOfTwo(nfft) then nfft / 2 else 0
        fullWithKernelSpectrumArrays(
          signal,
          kr,
          ki,
          m,
          blockLength,
          new Array[Double](nfft),
          new Array[Double](nfft),
          new Array[Double](packN),
          new Array[Double](packN),
          new Array[Double](nfft)
        )
      }

  /** Planned overlap-add with caller-owned FFT and block-output buffers. */
  def fullWithKernelSpectrumArrays(
      signal: DVec,
      kernelHalfRe: Array[Double],
      kernelHalfIm: Array[Double],
      kernelLength: Int,
      blockLength: Int,
      workRe: Array[Double],
      workIm: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double],
      blockOut: Array[Double]
  ): Either[SignalError, DVec] =
    val out = new Array[Double](signal.length + kernelLength - 1)
    fullWithKernelSpectrumArraysInto(
      signal,
      kernelHalfRe,
      kernelHalfIm,
      kernelLength,
      blockLength,
      workRe,
      workIm,
      packRe,
      packIm,
      blockOut,
      out
    ).map(_ => FftConvolution.adoptFull(out, out.length))

  /** Accumulates the full overlap-add result into caller-owned storage. */
  private[fft] def fullWithKernelSpectrumArraysInto(
      signal: DVec,
      kernelHalfRe: Array[Double],
      kernelHalfIm: Array[Double],
      kernelLength: Int,
      blockLength: Int,
      workRe: Array[Double],
      workIm: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double],
      blockOut: Array[Double],
      out: Array[Double]
  ): Either[SignalError, Unit] =
    fullWithKernelSpectrumArraysIntoRegion(
      signal,
      kernelHalfRe,
      kernelHalfIm,
      kernelLength,
      blockLength,
      fullOutLen = signal.length + kernelLength - 1,
      regionOffset = 0,
      regionLength = signal.length + kernelLength - 1,
      workRe,
      workIm,
      packRe,
      packIm,
      blockOut,
      out
    )

  /** Accumulates a requested full-overlap-add region into caller-owned storage. */
  private[fft] def fullWithKernelSpectrumArraysIntoRegion(
      signal: DVec,
      kernelHalfRe: Array[Double],
      kernelHalfIm: Array[Double],
      kernelLength: Int,
      blockLength: Int,
      fullOutLen: Int,
      regionOffset: Int,
      regionLength: Int,
      workRe: Array[Double],
      workIm: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double],
      blockOut: Array[Double],
      out: Array[Double]
  ): Either[SignalError, Unit] =
    if blockLength <= 0 then Left(SignalError.InvalidInputLength(blockLength))
    else
      val nfft = (kernelHalfRe.length - 1) << 1
      if blockLength + kernelLength - 1 > nfft then
        Left(SignalError.LengthMismatch(nfft, blockLength + kernelLength - 1))
      else if blockOut.length < nfft then
        Left(SignalError.LengthMismatch(nfft, blockOut.length))
      else if out.length < regionLength then
        Left(SignalError.LengthMismatch(regionLength, out.length))
      else if
        fullOutLen != signal.length + kernelLength - 1 ||
          regionOffset < 0 ||
          regionLength < 0 ||
          regionOffset + regionLength > fullOutLen
      then
        Left(
          SignalError.NumericalFailure(
            "OverlapAddConvolution",
            s"invalid output region offset=$regionOffset length=$regionLength fullLength=$fullOutLen"
          )
        )
      else
        java.util.Arrays.fill(out, 0, regionLength, 0.0)
        var pos = 0
        var err: Option[SignalError] = None
        while pos < signal.length && err.isEmpty do
          val len = math.min(blockLength, signal.length - pos)
          val blockOutLen = len + kernelLength - 1
          FftConvolution.fullBlockWithKernelSpectrumArraysInto(
            signal,
            pos,
            len,
            kernelHalfRe,
            kernelHalfIm,
            blockOutLen,
            workRe,
            workIm,
            packRe,
            packIm,
            blockOut
          ) match
            case Left(e) => err = Some(e)
            case Right(_) =>
              var i = 0
              while i < blockOutLen do
                val dest = pos + i
                if dest >= regionOffset && dest < regionOffset + regionLength then
                  out(dest - regionOffset) += blockOut(i)
                i += 1
          pos += blockLength
        err match
          case Some(e) => Left(e)
          case None    => Right(())
