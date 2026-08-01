package signal4s

extension [A](either: Either[SignalError, A])
  /** Convenience for examples and tests. Prefer matching on [[Left]] in library code. */
  def orThrow: A =
    either match
      case Right(value) => value
      case Left(error)  => throw new IllegalArgumentException(error.toString)
