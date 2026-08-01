package signal4s.laws

import gale.linalg.{DVec, Vec}
import signal4s.*
import signal4s.filter.*

class FilterLawSuite extends munit.FunSuite:

  private val x = Vec.tabulate(32)(i => math.sin(0.3 * i) + 0.1 * i)
  private val y = Vec.tabulate(17)(i => math.cos(0.2 * i))

  test("FIR batch from rest equals Input(Zero) convolution"):
    val fir = Fir.causal(Vec(0.2, 0.5, 0.3)).orThrow
    FilterLaws.firBatchEqualsConvolution(fir, x)

  test("FIR chunk concatenation and snapshot restart"):
    val fir = Fir.causal(Vec(1.0, -0.5, 0.25, 0.1)).orThrow
    val zero = FilterState.zeros(fir.stateLength).orThrow
    def step(chunk: DVec, st: FilterState): (DVec, FilterState) =
      val runner = fir.newRunner(st).orThrow
      val out = runner.process(chunk).orThrow
      (out, runner.snapshot)
    FilterLaws.chunkConcatenation(step, zero, x, y)
    FilterLaws.snapshotRestart(step, zero, x, y)

  test("FIR runners from one description have independent state"):
    val fir = Fir.causal(Vec(0.5, 0.5)).orThrow
    FilterLaws.independentRunners(() => fir.newRunner(), x)

  test("SOS chunk concatenation and snapshot restart"):
    val bq = Biquad(0.2, 0.4, 0.2, -0.5, 0.25).orThrow
    val sos = SecondOrderCascade.of(bq, bq).orThrow
    val zero = FilterState.zeros(sos.stateLength).orThrow
    def step(chunk: DVec, st: FilterState): (DVec, FilterState) =
      val runner = sos.newRunner(st).orThrow
      val out = runner.process(chunk).orThrow
      (out, runner.snapshot)
    FilterLaws.chunkConcatenation(step, zero, x, y)
    FilterLaws.snapshotRestart(step, zero, x, y)

  test("TF chunk concatenation"):
    val tf =
      DigitalTransferFunction(Vec(0.1, 0.2, 0.1), Vec(1.0, -0.5, 0.2)).orThrow
    val zero = FilterState.zeros(tf.stateLength).orThrow
    def step(chunk: DVec, st: FilterState): (DVec, FilterState) =
      val runner = tf.newRunner(st).orThrow
      val out = runner.process(chunk).orThrow
      (out, runner.snapshot)
    FilterLaws.chunkConcatenation(step, zero, x, y)
