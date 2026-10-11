package signal4s.fft

import gale.linalg.DVec
import signal4s.*

class BoundedStftParitySuite extends munit.FunSuite:
  test("pinned public SciPy full-frame complex coefficients match every frozen case"):
    for c <- BoundedStftGoldens.cases do
      val window = Window.fromTaps(DVec.fromSeq(c.window), WindowConvention.Periodic).orThrow
      val plan = BoundedStftPlan(window, SampleRate.hertz(16.0).orThrow, c.nfft, c.detrend, c.scaling).orThrow
      val output = plan.analyze(DVec.fromSeq(c.input), plan.resources.peakBytes).orThrow
      val scale = (c.real ++ c.imaginary).map(math.abs).max
      val tolerance = 1e-10 * scale + 1e-300
      assertEquals(output.bins.length, c.real.length)
      for k <- c.real.indices do
        assert(math.abs(output.bins.real(k) - c.real(k)) <= tolerance, s"${c.id} real $k")
        assert(math.abs(output.bins.imaginary(k) - c.imaginary(k)) <= tolerance, s"${c.id} imaginary $k")
        assertEqualsDouble(output.frequencies.frequencyAt(k).hertz, k.toDouble * 16 / c.nfft, 1e-12)
