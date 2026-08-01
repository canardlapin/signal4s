package signal4s

import gale.linalg.DVec

/** Reusable convolution algorithm selection. Implementations live in `signal4s-fft`. */
trait ConvolutionPlan:
  def kernel: Kernel
  def inputLength: Int
  def region: OutputRegion
  def selectedMethod: ConvolutionMethod
  def apply(signal: DVec): Either[SignalError, DVec]
