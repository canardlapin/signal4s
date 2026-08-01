package signal4s

/** Sealed error algebra for smart constructors and predictable numerical
  * failures. Successful but qualified calculations attach diagnostics to their
  * result types rather than expanding this enum.
  */
enum SignalError:
  case InvalidSampleRate(value: Double)
  case InvalidFrequency(value: Double)
  case InvalidDuration(value: Double)
  case FrequencyAboveNyquist(frequency: Frequency, sampleRate: SampleRate)
  case EmptyKernel
  case InvalidKernelOrigin(origin: Int, length: Int)
  case KernelNotOddLength(length: Int)
  case EmptySignal
  case EmptyConvolution
  case LengthMismatch(expected: Int, actual: Int)
  case InvalidPeriod(period: Int)
  case InvalidInputLength(value: Int)
  case UnsupportedOperatorRegion(region: OutputRegion)
  case InvalidFilterCoefficient(name: String, value: Double)
  case InvalidFilterState(expected: Int, actual: Int)
  case NumericalFailure(operation: String, detail: String)

  override def toString: String =
    this match
      case InvalidSampleRate(value) =>
        s"invalid sample rate: $value (must be finite and positive)"
      case InvalidFrequency(value) =>
        s"invalid frequency: $value (must be finite and non-negative)"
      case InvalidDuration(value) =>
        s"invalid duration: $value (must be finite)"
      case FrequencyAboveNyquist(frequency, sampleRate) =>
        s"frequency ${frequency.hertz} Hz exceeds Nyquist for sample rate ${sampleRate.hertz} Hz"
      case EmptyKernel =>
        "kernel must contain at least one tap"
      case InvalidKernelOrigin(origin, length) =>
        s"kernel origin $origin is outside tap range [0, $length)"
      case KernelNotOddLength(length) =>
        s"centeredOdd requires an odd tap count, got $length"
      case EmptySignal =>
        "signal must contain at least one sample"
      case EmptyConvolution =>
        "convolution produced an empty result for the requested region"
      case LengthMismatch(expected, actual) =>
        s"length mismatch: expected $expected, got $actual"
      case InvalidPeriod(period) =>
        s"invalid circular period: $period (must be positive)"
      case InvalidInputLength(value) =>
        s"invalid input length: $value (must be positive)"
      case UnsupportedOperatorRegion(region) =>
        s"convolution operator adjoint is only implemented for zero-extension regions; got $region"
      case InvalidFilterCoefficient(name, value) =>
        s"invalid filter coefficient $name=$value"
      case InvalidFilterState(expected, actual) =>
        s"filter state length mismatch: expected $expected, got $actual"
      case NumericalFailure(operation, detail) =>
        s"numerical failure in $operation: $detail"
