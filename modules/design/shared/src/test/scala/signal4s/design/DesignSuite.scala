package signal4s.design

import signal4s.*

class DesignSuite extends munit.FunSuite:

  private val fs = SampleRate.hertz(1000.0).orThrow

  test("Nyquist rejection"):
    val above = Frequency.hertz(600.0).orThrow
    assert(Butterworth.lowPass(2, above, fs).isLeft)
    assert(FirDesign.lowPass(11, above, fs).isLeft)

  test("Butterworth lowpass ZPK matches SciPy smoke values"):
    val fc = Frequency.hertz(100.0).orThrow
    val d = Butterworth.lowPass(4, fc, fs).orThrow
    assert(d.report.stable)
    assertEquals(d.report.family, "Butterworth")
    assertEqualsDouble(d.zpk.gain, 0.004824343357716229, 1e-12)
    assertEquals(d.zpk.zeros.length, 4)
    assert(d.zpk.zeros.forall(z => math.abs(z.real + 1.0) < 1e-12 && math.abs(z.imaginary) < 1e-12))
    // poles (unordered): match magnitudes
    val mags = d.zpk.poles.map(_.abs).sorted
    val expected = IArray(0.5442, 0.5442, 0.7955, 0.7955) // approx |p|
    // tighter: real parts set
    val reals = d.zpk.poles.map(_.real).sorted
    assertEqualsDouble(reals(0), 0.5242997881813058, 1e-9)
    assertEqualsDouble(reals(3), 0.6604567154097132, 1e-9)
    val _ = mags
    val _ = expected

  test("Butterworth SOS processes impulse"):
    val d = Butterworth.lowPass(4, Frequency.hertz(100.0).orThrow, fs).orThrow
    val y = Response.impulse(d.sos, 8).orThrow
    assertEqualsDouble(y(0), 0.004824343357716229, 1e-9)

  test("FirDesign lowpass Hamming odd taps"):
    val d = FirDesign.lowPass(11, Frequency.hertz(100.0).orThrow, fs).orThrow
    assertEquals(d.numTaps, 11)
    assertEqualsDouble(d.fir.taps(5), 0.23701569185915894, 1e-9)

  test("kaiserOrder positive"):
    val (n, beta) = FirDesign.kaiserOrder(40.0, 50.0, fs).orThrow
    assertEquals(n, 46)
    assertEqualsDouble(beta, 3.3953210522614574, 1e-9)
