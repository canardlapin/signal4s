package signal4s.laws

import gale.linalg.DVec
import signal4s.*
import signal4s.fft.*
import munit.Assertions

object FftLaws extends Assertions:

  def assertClose(a: ComplexVector, b: ComplexVector, tol: Double = 1e-9): Unit =
    assertEquals(a.length, b.length)
    var i = 0
    while i < a.length do
      val scale =
        math.max(
          1.0,
          math.max(
            math.hypot(a.real(i), a.imaginary(i)),
            math.hypot(b.real(i), b.imaginary(i))
          )
        )
      assert(
        math.hypot(a.real(i) - b.real(i), a.imaginary(i) - b.imaginary(i)) <= tol * scale,
        s"complex mismatch at $i"
      )
      i += 1

  def assertCloseReal(a: DVec, b: DVec, tol: Double = 1e-9): Unit =
    assertEquals(a.length, b.length)
    var i = 0
    while i < a.length do
      val scale = math.max(1.0, math.max(math.abs(a(i)), math.abs(b(i))))
      assert(math.abs(a(i) - b(i)) <= tol * scale, s"real mismatch at $i: ${a(i)} vs ${b(i)}")
      i += 1

  def complexRoundTrip(plan: FftPlan, x: ComplexVector): Unit =
    val y = plan.forward(x).orThrow
    val z = plan.inverse(y).orThrow
    assertClose(z, x, tol = 1e-8)

  def realRoundTrip(plan: RealFftPlan, x: DVec): Unit =
    val y = plan.forward(x).orThrow
    val z = plan.inverse(y).orThrow
    assertCloseReal(z, x, tol = 1e-8)

  /** Parseval: energy relation under the plan's normalization. */
  def parseval(plan: FftPlan, x: ComplexVector): Unit =
    val X = plan.forward(x).orThrow
    var timeE = 0.0
    var freqE = 0.0
    var i = 0
    while i < x.length do
      timeE += x.real(i) * x.real(i) + x.imaginary(i) * x.imaginary(i)
      freqE += X.real(i) * X.real(i) + X.imaginary(i) * X.imaginary(i)
      i += 1
    val expectedFreq =
      plan.normalization match
        case FftNormalization.Backward     => timeE * plan.length.toDouble
        case FftNormalization.Forward       => timeE / plan.length.toDouble
        case FftNormalization.Orthonormal   => timeE
    val scale = math.max(1.0, math.max(math.abs(expectedFreq), math.abs(freqE)))
    assert(math.abs(freqE - expectedFreq) <= 1e-7 * scale, s"Parseval $freqE vs $expectedFreq")

  def conjugateSymmetry(plan: RealFftPlan, x: DVec): Unit =
    val spec = plan.forward(x).orThrow
    val complexPlan = FftPlan(plan.length, plan.normalization).orThrow
    val cx = ComplexVector.tabulate(plan.length)(i => Complex(x(i), 0.0)).orThrow
    val fullSpec = complexPlan.forward(cx).orThrow
    var k = 1
    while k < plan.length do
      val a = fullSpec(k)
      val b = fullSpec(plan.length - k)
      assert(
        math.hypot(a.real - b.real, a.imaginary + b.imaginary) <= 1e-8 *
          math.max(1.0, a.abs),
        s"conjugate symmetry failed at k=$k"
      )
      k += 1
    var j = 0
    while j < spec.bins.length do
      assert(
        math.hypot(
          spec.bins.real(j) - fullSpec.real(j),
          spec.bins.imaginary(j) - fullSpec.imaginary(j)
        ) <= 1e-8,
        s"rfft bin $j mismatch"
      )
      j += 1
