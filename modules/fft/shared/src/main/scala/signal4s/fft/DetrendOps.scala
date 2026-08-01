package signal4s.fft

import gale.linalg.DVec

private[fft] object DetrendOps:

  def apply(frame: DVec, mode: Detrend): DVec =
    mode match
      case Detrend.None => frame
      case Detrend.Mean =>
        val mean = sum(frame) / frame.length.toDouble
        DVec.tabulate(frame.length)(i => frame(i) - mean)
      case Detrend.Linear =>
        // Least-squares line fit on sample indices 0..N-1
        val n = frame.length.toDouble
        var sumX = 0.0
        var sumY = 0.0
        var sumXX = 0.0
        var sumXY = 0.0
        var i = 0
        while i < frame.length do
          val x = i.toDouble
          val y = frame(i)
          sumX += x
          sumY += y
          sumXX += x * x
          sumXY += x * y
          i += 1
        val denom = n * sumXX - sumX * sumX
        if denom == 0.0 then apply(frame, Detrend.Mean)
        else
          val slope = (n * sumXY - sumX * sumY) / denom
          val intercept = (sumY - slope * sumX) / n
          DVec.tabulate(frame.length)(j => frame(j) - (intercept + slope * j.toDouble))

  private def sum(x: DVec): Double =
    var s = 0.0
    var i = 0
    while i < x.length do
      s += x(i)
      i += 1
    s
