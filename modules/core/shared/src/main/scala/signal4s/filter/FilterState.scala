package signal4s.filter

import gale.linalg.DVec
import signal4s.SignalError

/** Immutable snapshot of a runner's delay registers. */
final case class FilterState private (values: DVec):
  def length: Int = values.length

object FilterState:
  def zeros(length: Int): Either[SignalError, FilterState] =
    if length < 0 then Left(SignalError.InvalidFilterState(0, length))
    else Right(new FilterState(DVec.zeros(length)))

  def from(values: DVec): FilterState =
    new FilterState(values.copy)

  private[filter] def unsafe(values: DVec): FilterState =
    new FilterState(values)
