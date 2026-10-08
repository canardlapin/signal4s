package signal4s.fft.internal

/** Unscaled complex DFT (forward \(e^{-2\pi i}\), inverse \(e^{+2\pi i}\)).
  *
  * Strategy: iterative radix-2 when `n` is a power of two; the JVM delegates
  * 2/3/5-smooth complex lengths to JTransforms, while the portable engine uses
  * recursive Cooley–Tukey mixed-radix (2/3/4/5); Bluestein handles other
  * lengths. Mixed-radix reuses thread-local part buffers (amortized) and
  * specialized radix-2/4 combines.
  *
  * Power-of-two transforms reuse cached bit-reversal and twiddle tables.
  */
private[fft] object FftEngine:

  def forward(re: Array[Double], im: Array[Double]): Unit =
    transform(re, im, inverse = false)

  def inverse(re: Array[Double], im: Array[Double]): Unit =
    transform(re, im, inverse = true)

  private[fft] def supportsRealLength(n: Int): Boolean =
    n >= 2 && (isPowerOfTwo(n) || PlatformFft.supportsReal(n))

  private[fft] def hasNativeRealBackend: Boolean =
    PlatformFft.hasNativeRealBackend

  /** Unscaled real-input DFT via an \(n/2\) complex FFT (power-of-two `n ≥ 2`).
    *
    * On entry `re` holds the real signal; `im` is overwritten. On exit `re`/`im`
    * hold the full complex spectrum (Hermitian). `packRe`/`packIm` are length
    * `n/2` scratch buffers.
    */
  def forwardReal(
      re: Array[Double],
      im: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double]
  ): Unit =
    forwardRealCore(re, im, packRe, packIm, onesided = false)

  /** Like [[forwardReal]] but only fills bins `0 .. n/2` (Nyquist inclusive). */
  def forwardRealOnesided(
      re: Array[Double],
      im: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double]
  ): Unit =
    forwardRealCore(re, im, packRe, packIm, onesided = true)

  private def forwardRealCore(
      re: Array[Double],
      im: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double],
      onesided: Boolean
  ): Unit =
    val n = re.length
    require(n >= 2, "forwardReal requires length >= 2")
    require(im.length == n && packRe.length == n / 2 && packIm.length == n / 2)
    if PlatformFft.realForward(re, im, onesided) then return
    require(isPowerOfTwo(n), "portable forwardReal requires a power-of-two length")
    forwardRealPortable(re, im, packRe, packIm, onesided)

  private[fft] def forwardRealPortableForTest(
      re: Array[Double],
      im: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double],
      onesided: Boolean
  ): Unit =
    val n = re.length
    require(n >= 2 && isPowerOfTwo(n), "forwardReal requires power-of-two length >= 2")
    require(im.length == n && packRe.length == n / 2 && packIm.length == n / 2)
    forwardRealPortable(re, im, packRe, packIm, onesided)

  private def forwardRealPortable(
      re: Array[Double],
      im: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double],
      onesided: Boolean
  ): Unit =
    val n = re.length
    val N = n >> 1
    var k = 0
    while k < N do
      packRe(k) = re(2 * k)
      packIm(k) = re(2 * k + 1)
      k += 1
    radix2InPlace(packRe, packIm, inverse = false)
    val tw = RealUnpackTables.forLength(n)
    re(0) = packRe(0) + packIm(0)
    im(0) = 0.0
    re(N) = packRe(0) - packIm(0)
    im(N) = 0.0
    k = 1
    while k < N do
      val zr = packRe(k)
      val zi = packIm(k)
      val znr = packRe(N - k)
      val zni = packIm(N - k)
      val er = 0.5 * (zr + znr)
      val ei = 0.5 * (zi - zni)
      val or_ = 0.5 * (zr - znr)
      val oi = 0.5 * (zi + zni)
      val wr = tw.wRe(k)
      val wi = tw.wIm(k)
      val worR = wr * or_ - wi * oi
      val worI = wr * oi + wi * or_
      val xr = er + worI
      val xi = ei - worR
      re(k) = xr
      im(k) = xi
      if !onesided then
        re(n - k) = xr
        im(n - k) = -xi
      k += 1

  /** Unscaled inverse of [[forwardReal]]: Hermitian spectrum in → `n · x` in `re`. */
  def inverseReal(
      re: Array[Double],
      im: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double]
  ): Unit =
    if !PlatformFft.realInverse(re, im) then inverseRealPortable(re, im, packRe, packIm)

  private[fft] def inverseRealPortableForTest(
      re: Array[Double],
      im: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double]
  ): Unit =
    inverseRealPortable(re, im, packRe, packIm)

  private def inverseRealPortable(
      re: Array[Double],
      im: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double]
  ): Unit =
    val n = re.length
    inverseRealCore(re, im, packRe, packIm)
    val N = n >> 1
    var k = 0
    while k < N do
      re(2 * k) = 2.0 * packRe(k)
      re(2 * k + 1) = 2.0 * packIm(k)
      im(2 * k) = 0.0
      im(2 * k + 1) = 0.0
      k += 1

  /** Inverse real FFT writing `scale * (n · x)[0 until outLen]` into `out`.
    *
    * For convolution with Backward normalization, pass `scale = 1/n`.
    * Only bins `0 .. n/2` of `re`/`im` are read (Nyquist-inclusive onesided).
    */
  def inverseRealScaledTo(
      re: Array[Double],
      im: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double],
      scale: Double,
      out: Array[Double],
      outLen: Int
  ): Unit =
    val n = re.length
    require(outLen >= 0 && outLen <= n && out.length >= outLen)
    if PlatformFft.realInverse(re, im) then
      var i = 0
      while i < outLen do
        out(i) = scale * re(i)
        i += 1
    else
      inverseRealScaledToPortable(re, im, packRe, packIm, scale, out, outLen)

  private[fft] def inverseRealScaledToPortableForTest(
      re: Array[Double],
      im: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double],
      scale: Double,
      out: Array[Double],
      outLen: Int
  ): Unit =
    val n = re.length
    require(outLen >= 0 && outLen <= n && out.length >= outLen)
    inverseRealScaledToPortable(re, im, packRe, packIm, scale, out, outLen)

  private def inverseRealScaledToPortable(
      re: Array[Double],
      im: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double],
      scale: Double,
      out: Array[Double],
      outLen: Int
  ): Unit =
    inverseRealCore(re, im, packRe, packIm)
    val s = 2.0 * scale
    var i = 0
    while i + 1 < outLen do
      out(i) = s * packRe(i >> 1)
      out(i + 1) = s * packIm(i >> 1)
      i += 2
    if i < outLen then out(i) = s * packRe(i >> 1)

  /** Pack onesided Hermitian bins `0 .. n/2` into half-length `pack*` and inverse-FFT. */
  private def inverseRealCore(
      re: Array[Double],
      im: Array[Double],
      packRe: Array[Double],
      packIm: Array[Double]
  ): Unit =
    val n = re.length
    require(n >= 2 && isPowerOfTwo(n), "inverseReal requires power-of-two length >= 2")
    require(im.length == n && packRe.length == n / 2 && packIm.length == n / 2)
    val N = n >> 1
    val tw = RealUnpackTables.forLength(n)
    packRe(0) = 0.5 * (re(0) + re(N))
    packIm(0) = 0.5 * (re(0) - re(N))
    var k = 1
    while k < N do
      val xkr = re(k)
      val xki = im(k)
      // Pair bin k with bin N-k (both onesided); upper half n/2+1..n-1 is unused.
      val xnkr = re(N - k)
      val xnki = im(N - k)
      val er = 0.5 * (xkr + xnkr)
      val ei = 0.5 * (xki - xnki)
      val dr = xnkr - xkr
      val di = -xnki - xki
      val wr = tw.wRe(k)
      val wi = tw.wIm(k)
      val numR = di
      val numI = -dr
      val invWR = numR * wr + numI * wi
      val invWI = numI * wr - numR * wi
      packRe(k) = er + 0.5 * invWR
      packIm(k) = ei + 0.5 * invWI
      k += 1
    radix2InPlace(packRe, packIm, inverse = true)

  private def transform(re: Array[Double], im: Array[Double], inverse: Boolean): Unit =
    val n = re.length
    require(im.length == n, "real/imag length mismatch")
    if n <= 1 then ()
    else if isPowerOfTwo(n) then radix2InPlace(re, im, inverse)
    else if isSmooth(n) then
      if PlatformFft.complex(re, im, inverse) then ()
      else mixedRadixRecursive(re, im, inverse)
    else bluestein(re, im, inverse)

  /** Retains direct coverage of the portable smooth-length implementation when
    * a JVM acceleration backend is present.
    */
  private[fft] def mixedRadixPortableForTest(
      re: Array[Double],
      im: Array[Double],
      inverse: Boolean
  ): Unit =
    require(re.length == im.length, "real/imag length mismatch")
    require(re.length > 1 && !isPowerOfTwo(re.length) && isSmooth(re.length))
    mixedRadixRecursive(re, im, inverse)

  private def isPowerOfTwo(n: Int): Boolean =
    n > 0 && (n & (n - 1)) == 0

  private def isSmooth(n: Int): Boolean =
    var x = n
    while x % 2 == 0 do x /= 2
    while x % 3 == 0 do x /= 3
    while x % 5 == 0 do x /= 5
    x == 1

  /** Power-of-two FFT. Prefer a platform backend (JVM: JTransforms) when present;
    * otherwise bit-rev Cooley–Tukey. Stockham remains available for tests — it is
    * correct but was slower than bit-rev on Apple Silicon (2× working set).
    */
  private def radix2InPlace(re: Array[Double], im: Array[Double], inverse: Boolean): Unit =
    if PlatformFft.pow2(re, im, inverse) then ()
    else radix2Bitrev(re, im, inverse)

  /** Exposed for parity tests: Stockham vs bit-rev must match. */
  private[fft] def radix2StockhamForTest(
      re: Array[Double],
      im: Array[Double],
      tmpRe: Array[Double],
      tmpIm: Array[Double],
      inverse: Boolean
  ): Unit =
    radix2Stockham(re, im, tmpRe, tmpIm, inverse)

  private[fft] def radix2BitrevForTest(
      re: Array[Double],
      im: Array[Double],
      inverse: Boolean
  ): Unit =
    radix2Bitrev(re, im, inverse)

  /** In-place iterative radix-2 Cooley–Tukey with bit-reversal. */
  private def radix2Bitrev(re: Array[Double], im: Array[Double], inverse: Boolean): Unit =
    radix2Owned(re, im, inverse, Radix2Tables.forLength(re.length))

  /** Explicit tables: no platform scratch or global cache is touched. */
  private[fft] def radix2Owned(re: Array[Double], im: Array[Double], inverse: Boolean, tables: Radix2Tables): Unit =
    val n = re.length
    applyBitrev(re, im, tables.bitrev)
    val stages = if inverse then tables.invStages else tables.fwdStages
    var s = 0
    while s < stages.length do
      val stage = stages(s)
      val half = stage.wRe.length
      val len = half << 1
      val wReA = stage.wRe
      val wImA = stage.wIm
      var i = 0
      while i < n do
        var j = 0
        val half2 = half & ~1
        while j < half2 do
          val j0 = i + j
          val j1 = j0 + 1
          val h0 = j0 + half
          val h1 = j1 + half
          val w0r = wReA(j)
          val w0i = wImA(j)
          val w1r = wReA(j + 1)
          val w1i = wImA(j + 1)
          val u0r = re(j0)
          val u0i = im(j0)
          val u1r = re(j1)
          val u1i = im(j1)
          val v0r0 = re(h0)
          val v0i0 = im(h0)
          val v1r0 = re(h1)
          val v1i0 = im(h1)
          val v0r = v0r0 * w0r - v0i0 * w0i
          val v0i = v0r0 * w0i + v0i0 * w0r
          val v1r = v1r0 * w1r - v1i0 * w1i
          val v1i = v1r0 * w1i + v1i0 * w1r
          re(j0) = u0r + v0r
          im(j0) = u0i + v0i
          re(h0) = u0r - v0r
          im(h0) = u0i - v0i
          re(j1) = u1r + v1r
          im(j1) = u1i + v1i
          re(h1) = u1r - v1r
          im(h1) = u1i - v1i
          j += 2
        while j < half do
          val wRe = wReA(j)
          val wIm = wImA(j)
          val uRe = re(i + j)
          val uIm = im(i + j)
          val vRe0 = re(i + j + half)
          val vIm0 = im(i + j + half)
          val vRe = vRe0 * wRe - vIm0 * wIm
          val vIm = vRe0 * wIm + vIm0 * wRe
          re(i + j) = uRe + vRe
          im(i + j) = uIm + vIm
          re(i + j + half) = uRe - vRe
          im(i + j + half) = uIm - vIm
          j += 1
        i += len
      s += 1

  /** Stockham autosort radix-2 (Van Loan / scientificgo); result in `re`/`im`.
    *
    * Ping-pong buffers, no bit-reversal. Forward matches Go `s=+1`
    * (`W = e^{-iπ/l}`); inverse uses the conjugate recurrence.
    */
  private def radix2Stockham(
      re: Array[Double],
      im: Array[Double],
      tmpRe: Array[Double],
      tmpIm: Array[Double],
      inverse: Boolean
  ): Unit =
    val n = re.length
    require(tmpRe.length == n && tmpIm.length == n)
    val n2 = n >> 1
    // Match Go: y starts as a copy of x; each iteration swaps then writes y from tmp.
    System.arraycopy(re, 0, tmpRe, 0, n)
    System.arraycopy(im, 0, tmpIm, 0, n)
    var yRe: Array[Double] = tmpRe
    var yIm: Array[Double] = tmpIm
    var tRe: Array[Double] = re
    var tIm: Array[Double] = im

    // Go s=+1 forward → Sincos(-π/l); inverse s=-1 → Sincos(+π/l)
    val sTw = if inverse then -1.0 else 1.0
    var r = n2
    var l = 1
    while r >= 1 do
      val prevYRe = yRe
      val prevYIm = yIm
      yRe = tRe
      yIm = tIm
      tRe = prevYRe
      tIm = prevYIm

      val wAng = -sTw * math.Pi / l.toDouble
      val wre = math.cos(wAng)
      val wim = math.sin(wAng)
      var j = 0
      var wjRe = 1.0
      var wjIm = 0.0
      while j < l do
        val jrs = j * (r << 1)
        var k = jrs
        var m = jrs >> 1
        val kEnd = jrs + r
        while k < kEnd do
          val tr0 = tRe(k + r)
          val ti0 = tIm(k + r)
          val tr = wjRe * tr0 - wjIm * ti0
          val ti = wjRe * ti0 + wjIm * tr0
          val ur = tRe(k)
          val ui = tIm(k)
          yRe(m) = ur + tr
          yIm(m) = ui + ti
          yRe(m + n2) = ur - tr
          yIm(m + n2) = ui - ti
          m += 1
          k += 1
        val nwRe = wjRe * wre - wjIm * wim
        val nwIm = wjRe * wim + wjIm * wre
        wjRe = nwRe
        wjIm = nwIm
        j += 1
      l <<= 1
      r >>= 1

    if yRe ne re then
      System.arraycopy(yRe, 0, re, 0, n)
      System.arraycopy(yIm, 0, im, 0, n)

  private def applyBitrev(re: Array[Double], im: Array[Double], bitrev: Array[Int]): Unit =
    var i = 0
    while i < bitrev.length do
      val j = bitrev(i)
      if i < j then
        val tr = re(i); re(i) = re(j); re(j) = tr
        val ti = im(i); im(i) = im(j); im(j) = ti
      i += 1

  /** Recursive Cooley–Tukey for 2/3/5-smooth lengths (pooled buffers, cached twiddles). */
  private def mixedRadixRecursive(
      re: Array[Double],
      im: Array[Double],
      inverse: Boolean
  ): Unit =
    val n = re.length
    if n <= 1 then return
    if isPowerOfTwo(n) then
      radix2InPlace(re, im, inverse)
      return
    val radix =
      if n % 4 == 0 then 4
      else if n % 2 == 0 then 2
      else if n % 3 == 0 then 3
      else 5
    val m = n / radix
    val tw = MixedTwiddles.forLength(n, inverse)
    val parts = MixedParts.borrow(radix, m)
    try
      var r = 0
      while r < radix do
        val pRe = parts.re(r)
        val pIm = parts.im(r)
        var k = 0
        while k < m do
          pRe(k) = re(k * radix + r)
          pIm(k) = im(k * radix + r)
          k += 1
        mixedRadixRecursive(pRe, pIm, inverse)
        r += 1
      radix match
        case 2 => combineRadix2(re, im, parts, m, tw)
        case 4 => combineRadix4(re, im, parts, m, tw, inverse)
        case _ => combineGeneral(re, im, parts, radix, m, n, tw)
    finally MixedParts.release(parts)

  private def combineRadix2(
      re: Array[Double],
      im: Array[Double],
      parts: MixedParts,
      m: Int,
      tw: MixedTwiddles
  ): Unit =
    val y0r = parts.re(0)
    val y0i = parts.im(0)
    val y1r = parts.re(1)
    val y1i = parts.im(1)
    var k = 0
    while k < m do
      val wr = tw.wRe(k)
      val wi = tw.wIm(k)
      val a0r = y0r(k)
      val a0i = y0i(k)
      val b0r = y1r(k)
      val b0i = y1i(k)
      val t1r = b0r * wr - b0i * wi
      val t1i = b0r * wi + b0i * wr
      re(k) = a0r + t1r
      im(k) = a0i + t1i
      re(k + m) = a0r - t1r
      im(k + m) = a0i - t1i
      k += 1

  private def combineRadix4(
      re: Array[Double],
      im: Array[Double],
      parts: MixedParts,
      m: Int,
      tw: MixedTwiddles,
      inverse: Boolean
  ): Unit =
    val y0r = parts.re(0)
    val y0i = parts.im(0)
    val y1r = parts.re(1)
    val y1i = parts.im(1)
    val y2r = parts.re(2)
    val y2i = parts.im(2)
    val y3r = parts.re(3)
    val y3i = parts.im(3)
    val m2 = m << 1
    val m3 = m2 + m
    var k = 0
    while k < m do
      val w1r = tw.wRe(k)
      val w1i = tw.wIm(k)
      val k2 = k << 1
      val w2r = tw.wRe(k2)
      val w2i = tw.wIm(k2)
      val k3 = k2 + k
      val w3r = tw.wRe(k3)
      val w3i = tw.wIm(k3)
      val a0r = y0r(k)
      val a0i = y0i(k)
      val b1r = y1r(k)
      val b1i = y1i(k)
      val b2r = y2r(k)
      val b2i = y2i(k)
      val b3r = y3r(k)
      val b3i = y3i(k)
      val a1r = b1r * w1r - b1i * w1i
      val a1i = b1r * w1i + b1i * w1r
      val a2r = b2r * w2r - b2i * w2i
      val a2i = b2r * w2i + b2i * w2r
      val a3r = b3r * w3r - b3i * w3i
      val a3i = b3r * w3i + b3i * w3r
      val t0r = a0r + a2r
      val t0i = a0i + a2i
      val t1r = a0r - a2r
      val t1i = a0i - a2i
      val t2r = a1r + a3r
      val t2i = a1i + a3i
      val t3r = a1r - a3r
      val t3i = a1i - a3i
      // Multiply (t3r,t3i) by W_n^m = ±i (inverse: +i, forward: -i).
      val z3r = if inverse then -t3i else t3i
      val z3i = if inverse then t3r else -t3r
      re(k) = t0r + t2r
      im(k) = t0i + t2i
      re(k + m) = t1r + z3r
      im(k + m) = t1i + z3i
      re(k + m2) = t0r - t2r
      im(k + m2) = t0i - t2i
      re(k + m3) = t1r - z3r
      im(k + m3) = t1i - z3i
      k += 1

  private def combineGeneral(
      re: Array[Double],
      im: Array[Double],
      parts: MixedParts,
      radix: Int,
      m: Int,
      n: Int,
      tw: MixedTwiddles
  ): Unit =
    var k = 0
    while k < m do
      var t = 0
      while t < radix do
        var sumRe = 0.0
        var sumIm = 0.0
        var r = 0
        val base = t * m + k
        while r < radix do
          val idx = (base * r) % n
          val wr = tw.wRe(idx)
          val wi = tw.wIm(idx)
          val pr = parts.re(r)(k)
          val pi = parts.im(r)(k)
          sumRe += pr * wr - pi * wi
          sumIm += pr * wi + pi * wr
          r += 1
        re(base) = sumRe
        im(base) = sumIm
        t += 1
      k += 1

  private def bluestein(re: Array[Double], im: Array[Double], inverse: Boolean): Unit =
    val n = re.length
    val sign = if inverse then 1.0 else -1.0
    var m = 1
    while m < 2 * n - 1 do m <<= 1

    val aRe = Array.ofDim[Double](m)
    val aIm = Array.ofDim[Double](m)
    val bRe = Array.ofDim[Double](m)
    val bIm = Array.ofDim[Double](m)

    var i = 0
    while i < n do
      val angle = sign * math.Pi * i.toDouble * i.toDouble / n.toDouble
      val cr = math.cos(angle)
      val ci = math.sin(angle)
      aRe(i) = re(i) * cr - im(i) * ci
      aIm(i) = re(i) * ci + im(i) * cr
      bRe(i) = cr
      bIm(i) = -ci
      if i > 0 then
        bRe(m - i) = cr
        bIm(m - i) = -ci
      i += 1

    radix2InPlace(aRe, aIm, inverse = false)
    radix2InPlace(bRe, bIm, inverse = false)

    i = 0
    while i < m do
      val ar = aRe(i); val ai = aIm(i)
      val br = bRe(i); val bi = bIm(i)
      aRe(i) = ar * br - ai * bi
      aIm(i) = ar * bi + ai * br
      i += 1

    radix2InPlace(aRe, aIm, inverse = true)
    val invM = 1.0 / m.toDouble
    i = 0
    while i < n do
      val yr = aRe(i) * invM
      val yi = aIm(i) * invM
      val angle = sign * math.Pi * i.toDouble * i.toDouble / n.toDouble
      val cr = math.cos(angle)
      val ci = math.sin(angle)
      re(i) = yr * cr - yi * ci
      im(i) = yr * ci + yi * cr
      i += 1

/** Cached bit-reversal permutation and per-stage twiddles for one FFT length. */
private[fft] final class Radix2Tables private (
    val bitrev: Array[Int],
    val fwdStages: Array[Radix2Stage],
    val invStages: Array[Radix2Stage]
)

private[fft] final class Radix2Stage(
    val wRe: Array[Double],
    val wIm: Array[Double]
)

private[fft] object Radix2Tables:
  private var cache: Map[Int, Radix2Tables] = Map.empty

  def forLength(n: Int): Radix2Tables =
    cache.get(n) match
      case Some(t) => t
      case None =>
        val built = build(n)
        this.synchronized {
          cache.get(n) match
            case Some(t) => t
            case None =>
              cache = cache.updated(n, built)
              built
        }

  def owned(n: Int): Radix2Tables = build(n)

  private def build(n: Int): Radix2Tables =
    val bitrev = Array.ofDim[Int](n)
    var i = 0
    var j = 0
    bitrev(0) = 0
    i = 1
    while i < n do
      var bit = n >> 1
      while j >= bit do
        j -= bit
        bit >>= 1
      j += bit
      bitrev(i) = j
      i += 1

    var nStages = 0
    var len = 2
    while len <= n do
      nStages += 1
      len <<= 1

    val fwd = Array.ofDim[Radix2Stage](nStages)
    val inv = Array.ofDim[Radix2Stage](nStages)
    var s = 0
    len = 2
    while len <= n do
      val half = len / 2
      val fwdRe = Array.ofDim[Double](half)
      val fwdIm = Array.ofDim[Double](half)
      val invRe = Array.ofDim[Double](half)
      val invIm = Array.ofDim[Double](half)
      val angFwd = -2.0 * math.Pi / len.toDouble
      val angInv = 2.0 * math.Pi / len.toDouble
      var jj = 0
      while jj < half do
        val a = angFwd * jj.toDouble
        fwdRe(jj) = math.cos(a)
        fwdIm(jj) = math.sin(a)
        val b = angInv * jj.toDouble
        invRe(jj) = math.cos(b)
        invIm(jj) = math.sin(b)
        jj += 1
      fwd(s) = Radix2Stage(fwdRe, fwdIm)
      inv(s) = Radix2Stage(invRe, invIm)
      s += 1
      len <<= 1

    new Radix2Tables(bitrev, fwd, inv)

/** Twiddles \(W^k = e^{\pm 2\pi i k / n}\) for mixed-radix combine. */
private[internal] final class MixedTwiddles private (
    val wRe: Array[Double],
    val wIm: Array[Double]
)

private[internal] object MixedTwiddles:
  private var cache: Map[(Int, Boolean), MixedTwiddles] = Map.empty

  def forLength(n: Int, inverse: Boolean): MixedTwiddles =
    val key = (n, inverse)
    cache.get(key) match
      case Some(t) => t
      case None =>
        val built = build(n, inverse)
        this.synchronized {
          cache.get(key) match
            case Some(t) => t
            case None =>
              cache = cache.updated(key, built)
              built
        }

  private def build(n: Int, inverse: Boolean): MixedTwiddles =
    val sign = if inverse then 1.0 else -1.0
    val wRe = Array.ofDim[Double](n)
    val wIm = Array.ofDim[Double](n)
    var k = 0
    while k < n do
      val a = sign * 2.0 * math.Pi * k.toDouble / n.toDouble
      wRe(k) = math.cos(a)
      wIm(k) = math.sin(a)
      k += 1
    new MixedTwiddles(wRe, wIm)

/** Thread-local pool of mixed-radix part buffers (amortized zero alloc after warm-up). */
private[internal] final class MixedParts private (
    val re: Array[Array[Double]],
    val im: Array[Array[Double]],
    val partLen: Int
):
  val radix: Int = re.length

private[internal] object MixedParts:
  private val tl =
    new ThreadLocal[java.util.HashMap[Integer, java.util.ArrayDeque[MixedParts]]]:
      override def initialValue() =
        new java.util.HashMap[Integer, java.util.ArrayDeque[MixedParts]]()

  private def key(radix: Int, m: Int): Integer =
    Integer.valueOf((radix << 28) ^ m)

  def borrow(radix: Int, m: Int): MixedParts =
    val map = tl.get()
    val k = key(radix, m)
    var q = map.get(k)
    if q == null then
      q = new java.util.ArrayDeque[MixedParts]()
      val _ = map.put(k, q)
    val hit = q.pollFirst()
    if hit != null then hit
    else
      val re = Array.fill(radix)(Array.ofDim[Double](m))
      val im = Array.fill(radix)(Array.ofDim[Double](m))
      new MixedParts(re, im, m)

  def release(parts: MixedParts): Unit =
    val map = tl.get()
    val k = key(parts.radix, parts.partLen)
    var q = map.get(k)
    if q == null then
      q = new java.util.ArrayDeque[MixedParts]()
      val _ = map.put(k, q)
    val _ = q.addFirst(parts)

/** Twiddles \(W^k = e^{-2\pi i k / n}\) for real-FFT pack/unpack (`k = 0 .. n/2-1`). */
private[internal] final class RealUnpackTables private (
    val wRe: Array[Double],
    val wIm: Array[Double]
)

private[internal] object RealUnpackTables:
  private var cache: Map[Int, RealUnpackTables] = Map.empty

  def forLength(n: Int): RealUnpackTables =
    cache.get(n) match
      case Some(t) => t
      case None =>
        val built = build(n)
        this.synchronized {
          cache.get(n) match
            case Some(t) => t
            case None =>
              cache = cache.updated(n, built)
              built
        }

  private def build(n: Int): RealUnpackTables =
    val N = n >> 1
    val wRe = Array.ofDim[Double](N)
    val wIm = Array.ofDim[Double](N)
    var k = 0
    while k < N do
      val a = -2.0 * math.Pi * k.toDouble / n.toDouble
      wRe(k) = math.cos(a)
      wIm(k) = math.sin(a)
      k += 1
    new RealUnpackTables(wRe, wIm)
