package signal4s.fft

import gale.linalg.Vec
import signal4s.*

class RealWorkspaceContractSuite extends munit.FunSuite:
  test("forward refusal of same-length incompatible scratch preserves existing buffers"):
    val complex = FftPlan(4).orThrow
    val workspace = complex.newWorkspace()
    complex.forwardInto(ComplexVector(Vec(1.0, 2.0, -1.0, 0.0), Vec.zeros(4)).orThrow, workspace).orThrow
    val before = complex.resultFrom(workspace).orThrow
    val real = RealFftPlan(4, FftNormalization.Backward).orThrow
    assert(real.forwardInto(Vec(9.0, 8.0, 7.0, 6.0), workspace).isLeft)
    assertEquals(complex.resultFrom(workspace).orThrow.real.toSeq, before.real.toSeq)
    assertEquals(complex.resultFrom(workspace).orThrow.imaginary.toSeq, before.imaginary.toSeq)
    assertEquals(real.forwardInto(Vec.zeros(3), workspace), Left(SignalError.LengthMismatch(4, 3)))
    val good = real.newWorkspace()
    real.forwardInto(Vec(1.0, 2.0, -1.0, 0.0), good).orThrow
    assertEquals(real.spectrumFrom(good).orThrow.bins.real.toSeq, before.real.slice(0, 3).toSeq)

  test("inverse refusal preserves same-length incompatible buffers and accepts compatible scratch"):
    val real = RealFftPlan(4, FftNormalization.Backward).orThrow
    val input = Vec(1.0, 2.0, -1.0, 0.0)
    val spectrum = real.forward(input).orThrow
    val complex = FftPlan(4).orThrow
    val workspace = complex.newWorkspace()
    complex.forwardInto(ComplexVector(Vec(9.0, 8.0, 7.0, 6.0), Vec.zeros(4)).orThrow, workspace).orThrow
    val before = complex.resultFrom(workspace).orThrow
    assert(real.inverseInto(spectrum, workspace).isLeft)
    assertEquals(complex.resultFrom(workspace).orThrow.real.toSeq, before.real.toSeq)
    assertEquals(complex.resultFrom(workspace).orThrow.imaginary.toSeq, before.imaginary.toSeq)
    real.inverseInto(spectrum, real.newWorkspace()).orThrow
    assertEquals(real.inverse(spectrum).orThrow.toSeq, input.toSeq)
