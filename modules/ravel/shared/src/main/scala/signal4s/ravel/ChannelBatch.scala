package signal4s.ravel

import gale.linalg.{DVec, DMat, DMatBuilder}
import signal4s.SignalError

/** Planar multi-channel batch: one owned [[DVec]] per channel.
  *
  * Core 1-D types never depend on this; it lives only in `signal4s-ravel`.
  */
final case class ChannelBatch private (channels: IndexedSeq[DVec]):
  def channelCount: Int = channels.length
  def length: Int = if channels.isEmpty then 0 else channels.head.length

  def channel(i: Int): DVec = channels(i)

object ChannelBatch:
  def apply(channels: IndexedSeq[DVec]): Either[SignalError, ChannelBatch] =
    if channels.isEmpty then Right(new ChannelBatch(IndexedSeq.empty))
    else
      val n = channels.head.length
      var i = 1
      while i < channels.length do
        if channels(i).length != n then
          return Left(SignalError.LengthMismatch(n, channels(i).length))
        i += 1
      Right(new ChannelBatch(channels))

  def of(channels: DVec*): Either[SignalError, ChannelBatch] =
    apply(channels.toIndexedSeq)

  /** Pack planar channels into a matrix (rows = channels by default). */
  def toMatrix(batch: ChannelBatch, axis: Axis = Axis.Rows): Either[SignalError, DMat] =
    if batch.channelCount == 0 then Right(DMat.zeros(0, 0))
    else
      axis match
        case Axis.Rows =>
          val m = DMatBuilder.zeros(batch.channelCount, batch.length)
          var c = 0
          while c < batch.channelCount do
            var t = 0
            while t < batch.length do
              m(c, t) = batch.channels(c)(t)
              t += 1
            c += 1
          Right(m.result())
        case Axis.Columns =>
          val m = DMatBuilder.zeros(batch.length, batch.channelCount)
          var c = 0
          while c < batch.channelCount do
            var t = 0
            while t < batch.length do
              m(t, c) = batch.channels(c)(t)
              t += 1
            c += 1
          Right(m.result())

  def fromMatrix(mat: DMat, axis: Axis): Either[SignalError, ChannelBatch] =
    axis match
      case Axis.Rows =>
        val ch = IndexedSeq.tabulate(mat.rows)(r => mat.row(r).copy)
        apply(ch)
      case Axis.Columns =>
        val ch = IndexedSeq.tabulate(mat.cols)(c => mat.col(c).copy)
        apply(ch)
