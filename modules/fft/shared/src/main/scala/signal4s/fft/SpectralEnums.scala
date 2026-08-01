package signal4s.fft

/** Detrending applied to each frame before the FFT. */
enum Detrend:
  case None
  case Mean
  case Linear

/** Onesided (rfft) vs twosided spectrum. */
enum SpectralSides:
  case Onesided
  case Twosided

/** SciPy `scaling=` for periodogram / Welch. */
enum SpectralScaling:
  /** Power spectral density (power / Hz). */
  case Density
  /** Power spectrum (power per bin). */
  case Spectrum

/** How overlapping Welch segments are combined. */
enum AverageMethod:
  case Mean
  case Median
