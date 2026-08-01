package signal4s

import gale.linalg.Vec

class ConvolutionOperatorSuite extends munit.FunSuite:

  test("operator Full matches Convolution.apply"):
    val k = Kernel.causal(Vec(1.0, 2.0, 3.0)).orThrow
    val x = Vec(0.5, -1.0, 2.0, 0.0, 1.5)
    val op = Convolution.operator(k, x.length, OutputRegion.Full).orThrow
    val fromOp = op(x)
    val fromDirect = Convolution(x, k, OutputRegion.Full).orThrow
    assertEquals(fromOp.length, fromDirect.length)
    var i = 0
    while i < fromOp.length do
      assertEqualsDouble(fromOp(i), fromDirect(i), 1e-15)
      i += 1

  test("operator rejects Reflect boundary"):
    val k = Kernel.causal(Vec(1.0, 1.0)).orThrow
    assertEquals(
      Convolution.operator(k, 8, OutputRegion.Input(Boundary.Reflect)),
      Left(SignalError.UnsupportedOperatorRegion(OutputRegion.Input(Boundary.Reflect)))
    )

  test("operator rejects non-positive input length"):
    val k = Kernel.causal(Vec(1.0)).orThrow
    assertEquals(
      Convolution.operator(k, 0, OutputRegion.Full),
      Left(SignalError.InvalidInputLength(0))
    )

  test("small dense Toeplitz agreement for Full"):
    // y = convolve(x, [a,b]) => Toeplitz / banded structure
    val a = 2.0
    val b = -1.0
    val k = Kernel.causal(Vec(a, b)).orThrow
    val n = 4
    val op = Convolution.operator(k, n, OutputRegion.Full).orThrow
    // Manual matrix (rows = 5):
    // [a 0 0 0]
    // [b a 0 0]
    // [0 b a 0]
    // [0 0 b a]
    // [0 0 0 b]
    val x = Vec(1.0, 2.0, 3.0, 4.0)
    val y = op(x)
    assertEqualsDouble(y(0), a * 1.0, 1e-15)
    assertEqualsDouble(y(1), b * 1.0 + a * 2.0, 1e-15)
    assertEqualsDouble(y(2), b * 2.0 + a * 3.0, 1e-15)
    assertEqualsDouble(y(3), b * 3.0 + a * 4.0, 1e-15)
    assertEqualsDouble(y(4), b * 4.0, 1e-15)
    val g = Vec(1.0, 1.0, 1.0, 1.0, 1.0)
    val atg = op.adjoint(g)
    // A* 1 = [a+b, a+b, a+b, a+b] for this band pattern on interior...
    // (A*g)[k] = a*g[k] + b*g[k+1]
    assertEqualsDouble(atg(0), a * g(0) + b * g(1), 1e-15)
    assertEqualsDouble(atg(1), a * g(1) + b * g(2), 1e-15)
    assertEqualsDouble(atg(2), a * g(2) + b * g(3), 1e-15)
    assertEqualsDouble(atg(3), a * g(3) + b * g(4), 1e-15)

  test("valid centered operator has the direct adjoint"):
    val kernel = Kernel.centeredOdd(Vec(0.25, 0.5, 0.25)).orThrow
    val x = Vec(0.5, -1.0, 2.0, 0.0, 1.5, -0.25)
    val op = Convolution.operator(kernel, x.length, OutputRegion.Valid).orThrow
    val y = Vec(1.0, -0.5, 0.25, 2.0)
    val left = op(x).dot(y)
    val right = x.dot(op.adjoint(y))
    assertEqualsDouble(left, right, 1e-15)
