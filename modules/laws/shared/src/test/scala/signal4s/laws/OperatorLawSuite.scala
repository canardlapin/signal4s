package signal4s.laws

import gale.linalg.{DVec, Vec}
import signal4s.*

class OperatorLawSuite extends munit.ScalaCheckSuite:

  private val kernel = Kernel.causal(Vec(0.5, 1.0, -0.25, 0.125)).orThrow
  private val centered = Kernel.centeredOdd(Vec(0.25, 0.5, 0.25)).orThrow

  private def randomVec(n: Int, seed: Int): DVec =
    val rnd = new scala.util.Random(seed)
    DVec.tabulate(n)(_ => rnd.nextGaussian())

  test("full forward matches direct + adjoint law"):
    val op = Convolution.operator(kernel, inputLength = 17, OutputRegion.Full).orThrow
    val x = randomVec(17, seed = 1)
    val y = randomVec(op.rows, seed = 2)
    OperatorLaws.forwardMatchesDirect(op, x)
    OperatorLaws.adjointIdentity(op, x, y)
    OperatorLaws.agreesWithDenseMatrix(op)

  test("valid forward matches direct + adjoint law"):
    val op = Convolution.operator(kernel, inputLength = 20, OutputRegion.Valid).orThrow
    val x = randomVec(20, seed = 3)
    val y = randomVec(op.rows, seed = 4)
    OperatorLaws.forwardMatchesDirect(op, x)
    OperatorLaws.adjointIdentity(op, x, y)
    OperatorLaws.agreesWithDenseMatrix(op)

  test("input zero with nonzero origin: forward + adjoint"):
    val op =
      Convolution
        .operator(centered, inputLength = 13, OutputRegion.Input(Boundary.Zero))
        .orThrow
    val x = randomVec(13, seed = 5)
    val y = randomVec(op.rows, seed = 6)
    OperatorLaws.forwardMatchesDirect(op, x)
    OperatorLaws.adjointIdentity(op, x, y)
    OperatorLaws.agreesWithDenseMatrix(op)

  test("pathological short and odd lengths"):
    val cases = List(
      (3, Kernel.causal(Vec(1.0, 2.0)).orThrow, OutputRegion.Full),
      (5, Kernel.causal(Vec(1.0)).orThrow, OutputRegion.Full),
      (8, Kernel.centeredOdd(Vec(1.0, 0.0, -1.0)).orThrow, OutputRegion.Valid),
      (9, kernel, OutputRegion.Input(Boundary.Zero))
    )
    cases.zipWithIndex.foreach { case ((n, k, region), idx) =>
      val op = Convolution.operator(k, n, region).orThrow
      val x = randomVec(n, seed = 10 + idx)
      val y = randomVec(op.rows, seed = 20 + idx)
      OperatorLaws.forwardMatchesDirect(op, x)
      OperatorLaws.adjointIdentity(op, x, y)
    }

  test("non-zero boundaries are rejected"):
    OperatorLaws.rejectsNonZeroBoundary(kernel, inputLength = 10)
