package signal4s.laws

import gale.linalg.{DVec, Vec}
import signal4s.*
import org.scalacheck.Gen
import org.scalacheck.Prop.forAll

class ConvolutionLawSuite extends munit.ScalaCheckSuite:

  private val signal = Vec(0.5, -1.0, 2.0, 0.25, 3.0)
  private val taps = Vec(1.0, -0.5, 0.25)
  private val kernel = Kernel.causal(taps).orThrow
  private val finiteVector: Gen[DVec] =
    Gen
      .choose(1, 16)
      .flatMap(length => Gen.listOfN(length, Gen.chooseNum(-2.0, 2.0)))
      .map(DVec.fromSeq)

  test("delta identity"):
    ConvolutionLaws.deltaIdentity(signal)

  test("linearity"):
    val other = Vec(1.0, 1.0, 0.0, -2.0, 0.5)
    ConvolutionLaws.linearity(signal, other, kernel, alpha = 0.3, beta = -1.5)

  test("full commutativity"):
    ConvolutionLaws.fullCommutativity(signal, taps)

  test("full length"):
    ConvolutionLaws.fullLength(signal, kernel)

  test("signal full axis uses kernel origin"):
    val fs = SampleRate.hertz(200.0).orThrow
    val s = Signal(signal, Sampling(fs, Seconds.of(0.1).orThrow)).orThrow
    val centered = Kernel.centeredOdd(Vec(0.2, 0.6, 0.2)).orThrow
    ConvolutionLaws.signalFullAxis(s, centered)

  property("full convolution commutes across generated finite sequences"):
    forAll(finiteVector, finiteVector) { (x, h) =>
      ConvolutionLaws.fullCommutativity(x, h)
      true
    }
