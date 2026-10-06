package signal4s.fft

import signal4s.SignalError

/** Named window family with explicit length and convention. */
enum WindowSpec:
  case Rectangular(n: Int, conv: WindowConvention)
  case Hann(n: Int, conv: WindowConvention)
  case Hamming(n: Int, conv: WindowConvention)
  case Blackman(n: Int, conv: WindowConvention)
  case Kaiser(n: Int, beta: Double, conv: WindowConvention)

  def length: Int =
    this match
      case Rectangular(n, _) => n
      case Hann(n, _)        => n
      case Hamming(n, _)     => n
      case Blackman(n, _)    => n
      case Kaiser(n, _, _)   => n

  def convention: WindowConvention =
    this match
      case Rectangular(_, c) => c
      case Hann(_, c)        => c
      case Hamming(_, c)     => c
      case Blackman(_, c)    => c
      case Kaiser(_, _, c)   => c

  def validate: Either[SignalError, Unit] =
    if length <= 0 then Left(SignalError.InvalidInputLength(length))
    else
      this match
        case Kaiser(_, beta, _) if !beta.isFinite || beta < 0.0 =>
          Left(SignalError.NumericalFailure("WindowSpec", s"Kaiser beta must be >= 0, got $beta"))
        case _ => Right(())
