package signal4s.laws

import gale.linalg.{DVec, Vec}
import signal4s.*
import munit.Assertions

object ConvolutionLaws extends Assertions:

  def assertClose(actual: DVec, expected: DVec, tol: Double = 1e-10): Unit =
    assertEquals(actual.length, expected.length, "length mismatch")
    var i = 0
    while i < actual.length do
      val scale = math.max(1.0, math.max(math.abs(actual(i)), math.abs(expected(i))))
      assert(
        math.abs(actual(i) - expected(i)) <= tol * scale,
        s"mismatch at $i: ${actual(i)} vs ${expected(i)}"
      )
      i += 1

  /** Convolution with a unit impulse at lag 0 is the identity on Full. */
  def deltaIdentity(signal: DVec): Unit =
    val kernel = Kernel.causal(Vec(1.0)).orThrow
    val y = Convolution(signal, kernel, OutputRegion.Full).orThrow
    assertClose(y, signal)

  /** Convolution is linear in the signal. */
  def linearity(x: DVec, y: DVec, kernel: Kernel, alpha: Double, beta: Double): Unit =
    require(x.length == y.length)
    val axBy = DVec.tabulate(x.length)(i => alpha * x(i) + beta * y(i))
    val left = Convolution(axBy, kernel, OutputRegion.Full).orThrow
    val cx = Convolution(x, kernel, OutputRegion.Full).orThrow
    val cy = Convolution(y, kernel, OutputRegion.Full).orThrow
    val right = DVec.tabulate(left.length)(i => alpha * cx(i) + beta * cy(i))
    assertClose(left, right, tol = 1e-9)

  /** Full convolution commutes: x*h == h*x as coefficient sequences when both
    * are treated as causal kernels / signals of the same kind.
    */
  def fullCommutativity(x: DVec, h: DVec): Unit =
    val kx = Kernel.causal(x).orThrow
    val kh = Kernel.causal(h).orThrow
    val xy = Convolution(x, kh, OutputRegion.Full).orThrow
    val yx = Convolution(h, kx, OutputRegion.Full).orThrow
    assertClose(xy, yx, tol = 1e-9)

  def fullLength(signal: DVec, kernel: Kernel): Unit =
    val y = Convolution(signal, kernel, OutputRegion.Full).orThrow
    assertEquals(y.length, signal.length + kernel.length - 1)

  def signalFullAxis(signal: Signal, kernel: Kernel): Unit =
    val y = Convolution(signal, kernel, OutputRegion.Full).orThrow
    val expectedStart =
      signal.start.value - kernel.zeroLagIndex.toDouble / signal.sampleRate.hertz
    assertEqualsDouble(y.start.value, expectedStart, 1e-12)
    assertEquals(y.sampleRate.hertz, signal.sampleRate.hertz)
