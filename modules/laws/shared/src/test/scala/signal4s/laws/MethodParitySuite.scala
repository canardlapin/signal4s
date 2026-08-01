package signal4s.laws

import gale.linalg.Vec
import signal4s.*
import signal4s.fft.FftBackend

class MethodParitySuite extends munit.FunSuite:

  FftBackend.ensureInstalled()

  private def signal(n: Int) =
    Vec.tabulate(n)(i => math.sin(0.25 * i) + 0.05 * i)

  private def kernel(m: Int, origin: Int = 0) =
    val taps = Vec.tabulate(m)(i => 0.3 * math.cos(0.4 * i) + 0.1)
    Kernel.at(taps, origin).orThrow

  test("Direct ≡ FFT for Full/Valid/Input(Zero)"):
    val x = signal(48)
    val k = kernel(7, origin = 3)
    MethodParityLaws.directEqualsFft(x, k, OutputRegion.Full)
    MethodParityLaws.directEqualsFft(x, k, OutputRegion.Valid)
    MethodParityLaws.directEqualsFft(x, k, OutputRegion.Input(Boundary.Zero))

  test("Direct ≡ OLA for Full with several block sizes"):
    val x = signal(100)
    val k = kernel(9)
    List(8, 16, 32).foreach { b =>
      MethodParityLaws.directEqualsOla(x, k, b, OutputRegion.Full)
      MethodParityLaws.directEqualsOla(x, k, b, OutputRegion.Valid)
      MethodParityLaws.directEqualsOla(x, k, b, OutputRegion.Input(Boundary.Zero))
    }

  test("ConvolutionPlan reuses transformed kernel"):
    MethodParityLaws.planReusesKernel(signal(64), kernel(11, origin = 5))

  test("ownership stable after plan reuse"):
    MethodParityLaws.ownershipStable(signal(32), kernel(5))

  test("circular convolution theorem: Direct ≡ FFT product"):
    val x = signal(16)
    val k = kernel(5, origin = 2)
    MethodParityLaws.circularConvolutionTheorem(x, k)

  test("Auto resolves and matches Direct"):
    val x = signal(40)
    val k = kernel(6)
    val plan = Convolution.plan(k, x.length, OutputRegion.Full, ConvolutionMethod.Auto).orThrow
    val yAuto = plan(x).orThrow
    val yDirect = Convolution(x, k, OutputRegion.Full, ConvolutionMethod.Direct).orThrow
    ConvolutionLaws.assertClose(yAuto, yDirect, 1e-9)

  test("FFT rejected for non-zero Input boundary"):
    val x = signal(20)
    val k = kernel(3)
    val err = Convolution(x, k, OutputRegion.Input(Boundary.Reflect), ConvolutionMethod.Fft)
    assert(err.isLeft)

  test("large-length Direct ≡ FFT stress"):
    val x = signal(2048)
    val k = kernel(97, origin = 48)
    MethodParityLaws.directEqualsFft(x, k, OutputRegion.Full, tol = 1e-8)
    MethodParityLaws.directEqualsOla(x, k, 128, OutputRegion.Full, tol = 1e-8)
