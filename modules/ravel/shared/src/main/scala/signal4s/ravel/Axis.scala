package signal4s.ravel

/** Which matrix axis carries the 1-D signal samples. */
enum Axis:
  /** Each row is one signal (length = columns). */
  case Rows
  /** Each column is one signal (length = rows). */
  case Columns
