package signal4s.laws

import gale.linalg.{DVec, DoubleLinearOperator}
import signal4s.*
import munit.Assertions

object OperatorLaws extends Assertions:

  /** ⟨Ax, y⟩ ≈ ⟨x, A*y⟩ with a scale-aware tolerance. */
  def adjointIdentity(
      op: DoubleLinearOperator,
      x: DVec,
      y: DVec,
      tol: Double = 1e-9
  ): Unit =
    assertEquals(x.length, op.cols, "x length must equal op.cols")
    assertEquals(y.length, op.rows, "y length must equal op.rows")
    val ax = op(x)
    val aty = op.adjoint(y)
    val left = ax.dot(y)
    val right = x.dot(aty)
    val scale = math.max(1.0, math.max(math.abs(left), math.abs(right)))
    assert(
      math.abs(left - right) <= tol * scale,
      s"adjoint law failed: ⟨Ax,y⟩=$left ⟨x,A*y⟩=$right (tol=${tol * scale})"
    )

  /** Forward operator matches direct [[Convolution]] for the same region. */
  def forwardMatchesDirect(
      op: ConvolutionOperator,
      x: DVec
  ): Unit =
    assertEquals(x.length, op.inputLength)
    val fromOp = op(x)
    val fromDirect = Convolution(x, op.kernel, op.region).orThrow
    ConvolutionLaws.assertClose(fromOp, fromDirect, tol = 1e-12)

  /** Dense matrix from unit pulses agrees with matvec and adjoint matvec. */
  def agreesWithDenseMatrix(op: ConvolutionOperator): Unit =
    val dense = Array.tabulate(op.rows, op.cols) { (i, k) =>
      val e = DVec.tabulate(op.cols)(j => if j == k then 1.0 else 0.0)
      op(e)(i)
    }
    val x = DVec.tabulate(op.cols)(i => (i + 1).toDouble * 0.37 - 0.5)
    val y = DVec.tabulate(op.rows)(i => (i + 3).toDouble * 0.11 - 0.25)

    val ax = op(x)
    var i = 0
    while i < op.rows do
      var acc = 0.0
      var k = 0
      while k < op.cols do
        acc += dense(i)(k) * x(k)
        k += 1
      assert(
        math.abs(ax(i) - acc) <= 1e-12 * math.max(1.0, math.abs(acc)),
        s"forward dense mismatch at $i: ${ax(i)} vs $acc"
      )
      i += 1

    val aty = op.adjoint(y)
    var k = 0
    while k < op.cols do
      var acc = 0.0
      i = 0
      while i < op.rows do
        acc += dense(i)(k) * y(i)
        i += 1
      assert(
        math.abs(aty(k) - acc) <= 1e-12 * math.max(1.0, math.abs(acc)),
        s"adjoint dense mismatch at $k: ${aty(k)} vs $acc"
      )
      k += 1

  def rejectsNonZeroBoundary(kernel: Kernel, inputLength: Int): Unit =
    val err =
      Convolution.operator(
        kernel,
        inputLength,
        OutputRegion.Input(Boundary.Reflect)
      )
    assert(err.isLeft, "Reflect boundary must be rejected for operator")
