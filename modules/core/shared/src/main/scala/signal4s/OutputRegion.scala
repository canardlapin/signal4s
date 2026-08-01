package signal4s

/** Which output coordinates a convolution returns. */
enum OutputRegion:
  /** Entire finite-support convolution under zero extension. */
  case Full

  /** Positions for which every kernel tap overlaps observed input. */
  case Valid

  /** Exactly the input signal's coordinates, with the given boundary policy. */
  case Input(boundary: Boundary)
