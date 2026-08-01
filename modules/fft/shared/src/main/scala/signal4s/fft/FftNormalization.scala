package signal4s.fft

/** DFT scaling convention. The unscaled transform uses
  * \( X[k] = \sum_n x[n] e^{-2\pi i kn/N} \) (forward) and the conjugate
  * exponential without \(1/N\) (inverse).
  */
enum FftNormalization:
  /** Forward unscaled, inverse scaled by \(1/N\) (NumPy/SciPy default). */
  case Backward

  /** Forward scaled by \(1/N\), inverse unscaled. */
  case Forward

  /** Both directions scaled by \(1/\sqrt{N}\). */
  case Orthonormal
