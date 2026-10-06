package signal4s.multirate

import gale.linalg.{DVec, DVecBuilder}
import signal4s.SignalError

/** Upsample → FIR → downsample (SciPy `signal.upfirdn` semantics).
  *
  * Equivalent to zero-insertion by `up`, full convolution with `h`, then
  * keeping every `down`-th sample starting at index 0. Factors are GCD-reduced
  * through [[RateRatio]]; coefficients describe that reduced interpolation grid.
  * SciPy parity therefore requires the same reduced factors and prototype.
  */
object Upfirdn:

  def apply(
      h: DVec,
      x: DVec,
      up: Int,
      down: Int
  ): Either[SignalError, DVec] =
    RateRatio(up, down).flatMap(r => apply(h, x, r))

  def apply(h: DVec, x: DVec, ratio: RateRatio): Either[SignalError, DVec] =
    if h.length == 0 then Left(SignalError.EmptyKernel)
    else if x.length == 0 then Right(DVec.zeros(0))
    else if computedOutputLength(h.length, x.length, ratio.up, ratio.down) > Int.MaxValue then
      Left(SignalError.NumericalFailure("Upfirdn", "output length exceeds Int capacity"))
    else Right(polyphase(h, x, ratio.up, ratio.down))

  /** Output length for the supplied interpolation grid, matching SciPy `_output_len`.
    * To predict [[apply]], pass the GCD-reduced [[RateRatio]] factors. Invalid
    * dimensions or an unrepresentable Int result throw `IllegalArgumentException`.
    */
  def outputLength(hLen: Int, xLen: Int, up: Int, down: Int): Int =
    require(hLen > 0 && xLen >= 0 && up > 0 && down > 0,
      "output length requires positive taps/factors and nonnegative input length")
    val length = computedOutputLength(hLen, xLen, up, down)
    require(length <= Int.MaxValue, "output length exceeds Int capacity")
    length.toInt

  private def computedOutputLength(hLen: Int, xLen: Int, up: Int, down: Int): Long =
    if xLen == 0 then 0L
    else ((xLen - 1L) * up + hLen - 1L) / down + 1L

  /** Scatter only kept phases onto an Array; values match upsample→FIR→decimate. */
  private def polyphase(h: DVec, x: DVec, up: Int, down: Int): DVec =
    val outLen = outputLength(h.length, x.length, up, down)
    if outLen == 0 then return DVec.zeros(0)
    val m = h.length
    val n = x.length
    val hh = new Array[Double](m)
    val xx = new Array[Double](n)
    h.copyTo(hh)
    x.copyTo(xx)
    val acc = new Array[Double](outLen)

    if up == 1 && down == 1 then
      scatterFull(xx, hh, acc, n, m)
    else if down == 1 then
      scatterUpsampleOnly(xx, hh, acc, n, m, up, outLen)
    else if up == 1 && down <= math.min(m, n) then
      scatterDecimateBanked(xx, hh, acc, n, m, down, outLen)
    else
      scatterPolyphase(xx, hh, acc, n, m, up, down, outLen)

    val out = DVecBuilder.zeros(outLen)
    var i = 0
    while i < outLen do
      out(i) = acc(i)
      i += 1
    out.result()

  private def scatterFull(
      xx: Array[Double],
      hh: Array[Double],
      acc: Array[Double],
      n: Int,
      m: Int
  ): Unit =
    val m4 = m & ~3
    var xi = 0
    while xi < n do
      val xn = xx(xi)
      if xn != 0.0 then
        var k = 0
        while k < m4 do
          val base = xi + k
          acc(base) += xn * hh(k)
          acc(base + 1) += xn * hh(k + 1)
          acc(base + 2) += xn * hh(k + 2)
          acc(base + 3) += xn * hh(k + 3)
          k += 4
        while k < m do
          acc(xi + k) += xn * hh(k)
          k += 1
      xi += 1

  private def scatterUpsampleOnly(
      xx: Array[Double],
      hh: Array[Double],
      acc: Array[Double],
      n: Int,
      m: Int,
      up: Int,
      outLen: Int
  ): Unit =
    var xi = 0
    while xi < n do
      val xn = xx(xi)
      if xn != 0.0 then
        val base = xi * up
        val kMax = math.min(m, outLen - base)
        var k = 0
        val k4 = kMax & ~3
        while k < k4 do
          acc(base + k) += xn * hh(k)
          acc(base + k + 1) += xn * hh(k + 1)
          acc(base + k + 2) += xn * hh(k + 2)
          acc(base + k + 3) += xn * hh(k + 3)
          k += 4
        while k < kMax do
          acc(base + k) += xn * hh(k)
          k += 1
      xi += 1

  /** Precompute tap index lists per `xi % down` residue. */
  private def decimateTapBanks(m: Int, down: Int): Array[Array[Int]] =
    Array.tabulate(down) { r =>
      val start = (down - r) % down
      val nTaps = if start >= m then 0 else (m - 1 - start) / down + 1
      val ks = new Array[Int](nTaps)
      var t = 0
      var k = start
      while t < nTaps do
        ks(t) = k
        t += 1
        k += down
      ks
    }

  private def scatterDecimateBanked(
      xx: Array[Double],
      hh: Array[Double],
      acc: Array[Double],
      n: Int,
      m: Int,
      down: Int,
      outLen: Int
  ): Unit =
    val banks = decimateTapBanks(m, down)
    var xi = 0
    while xi < n do
      val xn = xx(xi)
      if xn != 0.0 then
        val ks = banks(xi % down)
        var t = 0
        while t < ks.length do
          val k = ks(t)
          val oi = (xi + k) / down
          if oi < outLen then acc(oi) += xn * hh(k)
          t += 1
      xi += 1

  private def scatterPolyphase(
      xx: Array[Double],
      hh: Array[Double],
      acc: Array[Double],
      n: Int,
      m: Int,
      up: Int,
      down: Int,
      outLen: Int
  ): Unit =
    val upQuot = up / down
    val upRem = up % down
    var xi = 0
    // For base = xi * up, maintain base = outBase * down + phase. This
    // removes the modulo at every input and the division at every tap.
    var outBase = 0
    var phase = 0
    while xi < n do
      val xn = xx(xi)
      if xn != 0.0 then
        var k = if phase == 0 then 0L else (down - phase).toLong
        var oi = outBase + (if phase == 0 then 0 else 1)
        while k < m && oi < outLen do
          acc(oi) += xn * hh(k.toInt)
          k += down.toLong
          oi += 1
      outBase += upQuot
      // Avoid overflowing Int when two individually valid residues are added.
      if phase >= down - upRem then
        phase -= down - upRem
        outBase += 1
      else phase += upRem
      xi += 1
