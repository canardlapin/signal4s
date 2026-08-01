package signal4s

/** Algorithm selection for convolution. Meaning is independent of method. */
enum ConvolutionMethod:
  case Auto
  case Direct
  case Fft
  case OverlapAdd(blockLength: Int)
