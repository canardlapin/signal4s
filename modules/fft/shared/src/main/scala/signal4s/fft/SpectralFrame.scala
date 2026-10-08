package signal4s.fft

import gale.linalg.DVec
import gale.numeric.ExactSum
import signal4s.*

/** Shared normalized frame arithmetic; the scale is restored after spectral normalization. */
private[fft] object SpectralFrame:
  def prepare(input: DVec, window: Window, detrend: Detrend): Either[SignalError, (DVec, Double)] =
    if input.length != window.length then Left(SignalError.LengthMismatch(window.length, input.length))
    else if !input.toSeq.forall(_.isFinite) then
      Left(SignalError.NumericalFailure("spectral.frame", "nonfinite input"))
    else
      val residual = detrend match
        case Detrend.None => Right((input, 1.0))
        case Detrend.Mean =>
          val anchor = input(0)
          val differences = DVec.tabulate(input.length)(i => input(i) - anchor)
          if !differences.toSeq.forall(_.isFinite) then
            return Left(SignalError.NumericalFailure("spectral.mean", "anchored difference exceeds finite capacity"))
          val total = ExactSum.zero()
          var i = 0
          while i < input.length do
            total.add(differences(i)) match
              case Left(e) => return Left(SignalError.NumericalFailure("spectral.mean", e.message))
              case Right(_) => ()
            i += 1
          val count = ExactSum.zero()
          count.add(input.length.toDouble) match
            case Left(e) => return Left(SignalError.NumericalFailure("spectral.mean", e.message))
            case Right(_) => ()
          total.ratio(count).left.map(e => SignalError.NumericalFailure("spectral.mean", e.message))
            .map(mean => (DVec.tabulate(input.length)(j => differences(j) - mean), 1.0))
        case Detrend.Linear =>
          val scale = input.toSeq.map(math.abs).max
          val normalized = if scale == 0 then input else DVec.tabulate(input.length)(j => input(j) / scale)
          Right((DetrendOps(normalized, detrend), scale))
      residual.flatMap { case (values, externalScale) =>
        if !values.toSeq.forall(_.isFinite) then
          Left(SignalError.NumericalFailure("spectral.detrend", "nonfinite detrended value"))
        else
          val magnitude = values.toSeq.map(math.abs).max
          val scale = magnitude * externalScale
          if !scale.isFinite then Left(SignalError.NumericalFailure("spectral.scale", "detrended magnitude exceeds finite capacity"))
          else
            val frame = DVec.tabulate(values.length)(i =>
              (if magnitude == 0 then 0.0 else values(i) / magnitude) * window.taps(i))
            if !frame.toSeq.forall(_.isFinite) then Left(SignalError.NumericalFailure("spectral.window", "nonfinite windowed value"))
            else Right((frame, scale))
      }

  /** Exponent arithmetic avoids overflow in an unnormalized square or scale product. */
  def power(re: Double, im: Double, scale: Double, divisor: Double, factor: Double): Either[SignalError, Double] =
    val magnitude = math.hypot(re, im)
    if !magnitude.isFinite then Left(SignalError.NumericalFailure("spectral.fft", "nonfinite FFT magnitude"))
    else if magnitude == 0 || scale == 0 then Right(0.0)
    else
      def parts(value: Double): (Double, Int) =
        val initial = java.lang.Math.getExponent(value)
        val exponent = if initial == -1023 then java.lang.Math.getExponent(java.lang.Math.scalb(value, 54)) - 54 else initial
        (java.lang.Math.scalb(value, -exponent), exponent)
      val (a, ae) = parts(magnitude)
      val (b, be) = parts(scale)
      val (c, ce) = parts(divisor)
      val value = java.lang.Math.scalb(a * a * b * b / c * factor, 2 * (ae + be) - ce)
      if value.isFinite then Right(value)
      else Left(SignalError.NumericalFailure("spectral.power", "normalized power exceeds finite Double capacity"))
