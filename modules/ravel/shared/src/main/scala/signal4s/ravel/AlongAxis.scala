package signal4s.ravel

import gale.linalg.{DVec, DMat, DMatBuilder}
import signal4s.SignalError

/** Apply a 1-D transform independently along matrix rows or columns. */
object AlongAxis:

  def map(
      mat: DMat,
      axis: Axis,
      policy: ContiguousPolicy = ContiguousPolicy.CopyAlways
  )(f: DVec => Either[SignalError, DVec]): Either[SignalError, DMat] =
    axis match
      case Axis.Rows =>
        mapRows(mat, policy)(f)
      case Axis.Columns =>
        mapColumns(mat, policy)(f)

  def mapChannels(
      batch: ChannelBatch,
      policy: ContiguousPolicy = ContiguousPolicy.CopyAlways
  )(f: DVec => Either[SignalError, DVec]): Either[SignalError, ChannelBatch] =
    if batch.channelCount == 0 then Right(batch)
    else
      val out = IndexedSeq.newBuilder[DVec]
      var i = 0
      var err: Option[SignalError] = None
      while i < batch.channelCount && err.isEmpty do
        val src =
          policy match
            case ContiguousPolicy.CopyAlways => batch.channel(i).copy
            case ContiguousPolicy.RequireContiguous =>
              // Planar channels are already owned contiguous vectors.
              batch.channel(i)
        f(src) match
          case Left(e)  => err = Some(e)
          case Right(y) => out += y
        i += 1
      err match
        case Some(e) => Left(e)
        case None    => ChannelBatch(out.result())

  private def mapRows(
      mat: DMat,
      policy: ContiguousPolicy
  )(f: DVec => Either[SignalError, DVec]): Either[SignalError, DMat] =
    if mat.rows == 0 then Right(mat)
    else
      val first = prepare(mat.row(0), mat.isContiguousRowMajor, policy).flatMap(f)
      first.flatMap { y0 =>
        val out = DMatBuilder.zeros(mat.rows, y0.length)
        var c = 0
        while c < y0.length do
          out(0, c) = y0(c)
          c += 1
        var r = 1
        var err: Option[SignalError] = None
        while r < mat.rows && err.isEmpty do
          prepare(mat.row(r), mat.isContiguousRowMajor, policy).flatMap(f) match
            case Left(e) => err = Some(e)
            case Right(y) =>
              if y.length != y0.length then
                err = Some(SignalError.LengthMismatch(y0.length, y.length))
              else
                var j = 0
                while j < y.length do
                  out(r, j) = y(j)
                  j += 1
          r += 1
        err match
          case Some(e) => Left(e)
          case None    => Right(out.result())
      }

  private def mapColumns(
      mat: DMat,
      policy: ContiguousPolicy
  )(f: DVec => Either[SignalError, DVec]): Either[SignalError, DMat] =
    if mat.cols == 0 then Right(mat)
    else
      val first = prepare(mat.col(0), mat.isContiguousColMajor, policy).flatMap(f)
      first.flatMap { y0 =>
        val out = DMatBuilder.zeros(y0.length, mat.cols)
        var r = 0
        while r < y0.length do
          out(r, 0) = y0(r)
          r += 1
        var c = 1
        var err: Option[SignalError] = None
        while c < mat.cols && err.isEmpty do
          prepare(mat.col(c), mat.isContiguousColMajor, policy).flatMap(f) match
            case Left(e) => err = Some(e)
            case Right(y) =>
              if y.length != y0.length then
                err = Some(SignalError.LengthMismatch(y0.length, y.length))
              else
                var i = 0
                while i < y.length do
                  out(i, c) = y(i)
                  i += 1
          c += 1
        err match
          case Some(e) => Left(e)
          case None    => Right(out.result())
      }

  private def prepare(
      view: DVec,
      contiguousAlongAxis: Boolean,
      policy: ContiguousPolicy
  ): Either[SignalError, DVec] =
    policy match
      case ContiguousPolicy.CopyAlways =>
        Right(view.copy)
      case ContiguousPolicy.RequireContiguous =>
        if contiguousAlongAxis then Right(view)
        else
          Left(
            SignalError.NumericalFailure(
              "AlongAxis",
              "RequireContiguous: axis storage is strided; use CopyAlways or gather first"
            )
          )
