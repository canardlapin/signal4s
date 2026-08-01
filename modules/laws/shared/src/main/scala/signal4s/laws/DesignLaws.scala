package signal4s.laws

import signal4s.*
import signal4s.design.*
import munit.Assertions

object DesignLaws extends Assertions:

  /** ZPK→TF and SOS frequency responses agree in magnitude. */
  def zpkSosTfAgree(designed: DesignedIir, sampleRate: SampleRate, nfft: Int = 64): Unit =
    val tf = Transfer.fromZerosPolesGain(designed.zpk).orThrow
    val hTf = Response.freqz(tf, sampleRate, nfft).orThrow.magnitude
    val hSos = Response.sosFreqz(designed.sos, sampleRate, nfft).orThrow.magnitude
    assertEquals(hTf.length, hSos.length)
    var i = 0
    while i < hTf.length do
      val scale = math.max(math.abs(hTf(i)), math.abs(hSos(i)))
      val tol = 1e-9 + 1e-6 * scale
      assert(
        math.abs(hTf(i) - hSos(i)) <= tol,
        s"mag mismatch at bin $i: ${hTf(i)} vs ${hSos(i)}"
      )
      i += 1

  def butterworthStable(order: Int, cutoffHz: Double, fsHz: Double): Unit =
    val fs = SampleRate.hertz(fsHz).orThrow
    val fc = Frequency.hertz(cutoffHz).orThrow
    val d = Butterworth.lowPass(order, fc, fs).orThrow
    assert(d.report.stable, s"unstable butterworth order=$order")
    assert(d.report.maxPoleMagnitude < 1.0)
