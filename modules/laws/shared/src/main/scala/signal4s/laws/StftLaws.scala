package signal4s.laws

import gale.linalg.DVec
import signal4s.*
import signal4s.fft.*
import munit.Assertions

object StftLaws extends Assertions:

  def roundTrip(plan: StftPlan, signal: DVec, tol: Double = 1e-9): Unit =
    val tf = plan.analyze(signal).orThrow
    val y = plan.synthesize(tf, signal.length).orThrow
    assertEquals(y.length, signal.length)
    var i = 0
    while i < signal.length do
      val scale = math.max(1.0, math.max(math.abs(signal(i)), math.abs(y(i))))
      assert(
        math.abs(y(i) - signal(i)) <= tol * scale,
        s"STFT round-trip mismatch at $i: ${y(i)} vs ${signal(i)}"
      )
      i += 1

  def axesMatchPlan(plan: StftPlan, signal: DVec): Unit =
    val tf = plan.analyze(signal).orThrow
    assertEquals(tf.nfft, plan.nfft)
    assertEquals(tf.times.length, tf.frameCount)
    assertEquals(tf.frequencies.length, plan.nfft / 2 + 1)
    assertEquals(tf.times.step.value, plan.hop.toDouble / plan.sampleRate.hertz, 1e-12)
    val _ = tf.spectrogram

  def dualComputedOnce(plan: StftPlan): Unit =
    assert(plan.dual.isDefined, "expected CanonicalDual at construction")
    assertEquals(plan.dual.get.length, plan.frameLength)
    assertEquals(plan.synthesis.length, plan.frameLength)
