package signal4s

/** How to read samples outside the observed support. */
enum Boundary:
  case Zero
  case Constant(value: Double)
  case Clamp
  /** Reflect without repeating the edge sample (NumPy `reflect`). */
  case Reflect
  /** Reflect repeating the edge sample (NumPy `symmetric`). */
  case Symmetric
