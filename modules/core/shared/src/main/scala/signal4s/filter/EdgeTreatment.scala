package signal4s.filter

/** Edge / initialization policy for offline forward–backward filtering. */
enum EdgeTreatment:
  /** Odd reflection about the endpoints (SciPy `filtfilt` `padtype='odd'`). */
  case OddPad

  /** Even reflection about the endpoints (SciPy `padtype='even'`). */
  case EvenPad

  /** Constant extension (SciPy `padtype='constant'`). */
  case ConstantPad(value: Double)

  /** Gustafsson's method (SciPy `filtfilt(..., method='gust')`). */
  case Gustafsson
