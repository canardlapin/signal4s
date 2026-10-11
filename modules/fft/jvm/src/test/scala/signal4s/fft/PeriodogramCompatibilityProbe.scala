package signal4s.fft

import gale.linalg.Vec
import signal4s.*
import java.nio.file.{Files, Path}

/** Same binary probe can run with the old locked or candidate FFT classes first. */
object PeriodogramCompatibilityProbe:
  def main(args: Array[String]): Unit =
    val output = Vector.newBuilder[String]
    val fs = SampleRate.hertz(16.0).orThrow
    for nfft <- Vector(5, 7, 8, 11, 16);
        detrend <- Vector(Detrend.None, Detrend.Mean, Detrend.Linear);
        scaling <- Vector(SpectralScaling.Density, SpectralScaling.Spectrum);
        hann <- Vector(false, true) do
      val window = Window.fromSpec(if hann then WindowSpec.Hann(5, WindowConvention.Periodic)
        else WindowSpec.Rectangular(5, WindowConvention.Periodic)).orThrow
      val input = Vec(1.0, -2.0, 0.5, 3.0, -0.25)
      val result = BoundedPeriodogramPlan(window, fs, nfft, detrend, scaling).orThrow.estimate(input).orThrow
      val payload = result.power.toSeq.map(v => java.lang.Long.toHexString(java.lang.Double.doubleToRawLongBits(v))).mkString(",")
      output += s"$nfft|$detrend|$scaling|$hann|$payload"
    Files.writeString(Path.of(args(0)), output.result().mkString("\n") + "\n")
    ()
