package signal4s.filter

import gale.linalg.{DVec, MutableDVec}

/** Direct-form II transposed kernels (SciPy `lfilter` / `sosfilt` convention). */
private[filter] object Df2Transposed:

  /** One biquad step. `state` length 2, mutated in place. */
  def biquadStep(section: Biquad, x: Double, state: Array[Double]): Double =
    val y = section.b0 * x + state(0)
    state(0) = section.b1 * x - section.a1 * y + state(1)
    state(1) = section.b2 * x - section.a2 * y
    y

  /** General TF step. `b`/`a` are delay-power coeffs with `a(0) == 1`.
    * `state` length is `max(b.length, a.length) - 1`.
    */
  def transferStep(
      b: DVec,
      a: DVec,
      x: Double,
      state: Array[Double]
  ): Double =
    val y = b(0) * x + (if state.length > 0 then state(0) else 0.0)
    val order = state.length
    var i = 0
    while i < order - 1 do
      val bi = if i + 1 < b.length then b(i + 1) else 0.0
      val ai = if i + 1 < a.length then a(i + 1) else 0.0
      state(i) = bi * x - ai * y + state(i + 1)
      i += 1
    if order > 0 then
      val bi = if order < b.length then b(order) else 0.0
      val ai = if order < a.length then a(order) else 0.0
      state(order - 1) = bi * x - ai * y
    y

  def filterInto(
      input: DVec,
      output: MutableDVec,
      step: Double => Double
  ): Unit =
    var n = 0
    while n < input.length do
      output(n) = step(input(n))
      n += 1
