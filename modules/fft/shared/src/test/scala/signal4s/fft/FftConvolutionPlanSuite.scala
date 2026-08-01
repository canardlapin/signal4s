package signal4s.fft

import gale.linalg.{DVec, Vec}
import signal4s.*

/** Exercise reusable FFT/OLA plan contracts as well as numerical happy paths. */
class FftConvolutionPlanSuite extends munit.FunSuite:

  FftBackend.ensureInstalled()

  private val signal = Vec(1.0, -2.0, 0.5, 3.0, -1.0)
  private val kernel = Kernel.centeredOdd(Vec(0.25, 0.5, 0.25)).orThrow

  test("kernel spectra reject invalid transform lengths and reuse correctly"):
    assertEquals(FftConvolution.transformKernel(kernel, 0), Left(SignalError.InvalidInputLength(0)))
    assertEquals(
      FftConvolution.transformKernel(kernel, 2),
      Left(SignalError.LengthMismatch(expected = 2, actual = 3))
    )

    val nfft = FftConvolution.fftLength(signal.length + kernel.length - 1)
    val spectrum = FftConvolution.transformKernel(kernel, nfft).orThrow
    val actual =
      FftConvolution.fullWithKernelSpectrum(signal, spectrum, signal.length + kernel.length - 1).orThrow
    val expected = Convolution(signal, kernel, OutputRegion.Full, ConvolutionMethod.Direct).orThrow
    assertClose(actual, expected)

  test("FFT work-buffer validation rejects every incompatible shape"):
    val halfRe = new Array[Double](5) // nfft = 8
    val halfIm = new Array[Double](5)
    val workRe = new Array[Double](8)
    val workIm = new Array[Double](8)
    val packRe = new Array[Double](4)
    val packIm = new Array[Double](4)
    val signal = Vec(1.0, 2.0)

    assertEquals(
      FftConvolution.fullWithKernelSpectrumArraysInto(
        signal,
        halfRe,
        new Array[Double](4),
        3,
        workRe,
        workIm,
        packRe,
        packIm
      ),
      Left(SignalError.LengthMismatch(expected = 5, actual = 4))
    )
    assertEquals(
      FftConvolution.fullWithKernelSpectrumArraysInto(
        signal,
        halfRe,
        halfIm,
        3,
        new Array[Double](7),
        workIm,
        packRe,
        packIm
      ),
      Left(SignalError.LengthMismatch(expected = 8, actual = 7))
    )
    assertEquals(
      FftConvolution.fullWithKernelSpectrumArraysInto(
        signal,
        halfRe,
        halfIm,
        3,
        workRe,
        workIm,
        packRe,
        packIm,
        new Array[Double](2)
      ),
      Left(SignalError.LengthMismatch(expected = 3, actual = 2))
    )
    assert(
      FftConvolution
        .fullWithKernelSpectrumArraysInto(
          signal,
          halfRe,
          halfIm,
          -1,
          workRe,
          workIm,
          packRe,
          packIm,
          new Array[Double](8)
        )
        .isLeft
    )
    assertEquals(
      FftConvolution.fullWithKernelSpectrumArraysInto(
        Vec.tabulate(9)(_.toDouble),
        halfRe,
        halfIm,
        3,
        workRe,
        workIm,
        packRe,
        packIm
      ),
      Left(SignalError.LengthMismatch(expected = 8, actual = 9))
    )
    assertEquals(
      FftConvolution.fullWithKernelSpectrumArraysInto(
        signal,
        halfRe,
        halfIm,
        3,
        workRe,
        workIm,
        new Array[Double](3),
        packIm
      ),
      Left(SignalError.LengthMismatch(expected = 4, actual = 3))
    )

  test("region extraction matches direct semantics and rejects nonzero boundaries"):
    val full = FftConvolution.full(signal, kernel).orThrow
    List(OutputRegion.Full, OutputRegion.Valid, OutputRegion.Input(Boundary.Zero)).foreach { region =>
      val actual = FftConvolution.extractRegion(full, signal.length, kernel, region).orThrow
      val expected = Convolution(signal, kernel, region, ConvolutionMethod.Direct).orThrow
      assertClose(actual, expected)
    }
    assert(
      FftConvolution
        .extractRegion(full, signal.length, kernel, OutputRegion.Input(Boundary.Reflect))
        .isLeft
    )

  test("planned FFT and overlap-add reject unsupported shapes and input lengths"):
    assert(
      Convolution.plan(
        kernel,
        signal.length,
        OutputRegion.Input(Boundary.Reflect),
        ConvolutionMethod.Fft
      ).isLeft
    )
    assert(OverlapAddConvolution.full(signal, kernel, blockLength = 0).isLeft)

    val plan = Convolution.plan(kernel, signal.length, OutputRegion.Input(Boundary.Zero), ConvolutionMethod.Fft).orThrow
    assertEquals(plan(Vec(1.0, 2.0)), Left(SignalError.LengthMismatch(signal.length, 2)))

    val ola = Convolution.plan(kernel, signal.length, OutputRegion.Valid, ConvolutionMethod.OverlapAdd(2)).orThrow
    val actual = ola(signal).orThrow
    val expected = Convolution(signal, kernel, OutputRegion.Valid, ConvolutionMethod.Direct).orThrow
    assertClose(actual, expected)

  private def assertClose(actual: DVec, expected: DVec): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assertEqualsDouble(actual(i), expected(i), 1e-10, clue = s"i=$i")
      i += 1
