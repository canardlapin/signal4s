package signal4s.fft.internal

import org.jtransforms.fft.DoubleFFT_1D
import pl.edu.icm.jlargearrays.ConcurrencyUtils
import java.util.concurrent.ConcurrentHashMap

/** JVM complex FFT via JTransforms (Ooura-derived, pure Java). */
private[internal] object PlatformFft:

  // Match signal4s single-owner plan discipline; SciPy pocketfft is typically 1-thread.
  ConcurrencyUtils.setNumberOfThreads(1)

  def hasNativeRealBackend: Boolean = true

  private val plans = new ConcurrentHashMap[Integer, DoubleFFT_1D]()

  private def plan(n: Int): DoubleFFT_1D =
    val key = Integer.valueOf(n)
    var p = plans.get(key)
    if p == null then
      p = new DoubleFFT_1D(n)
      val prev = plans.putIfAbsent(key, p)
      if prev != null then p = prev
    p

  /** @return true if handled. */
  def pow2(re: Array[Double], im: Array[Double], inverse: Boolean): Boolean =
    val n = re.length
    if n < 2 || (n & (n - 1)) != 0 || im.length != n then false
    else complex(re, im, inverse)

  /** JTransforms supports arbitrary complex lengths. The engine elects this
    * backend for JVM 2/3/5-smooth lengths and keeps portable Bluestein as the
    * fallback for other shapes.
    */
  def complex(re: Array[Double], im: Array[Double], inverse: Boolean): Boolean =
    val n = re.length
    if n < 2 || im.length != n then false
    else
      val buf = InterleavedScratch.borrow(n)
      try
        var i = 0
        while i < n do
          buf(2 * i) = re(i)
          buf(2 * i + 1) = im(i)
          i += 1
        val fft = plan(n)
        if inverse then fft.complexInverse(buf, false) // scale=false → unscaled like our engine
        else fft.complexForward(buf)
        i = 0
        while i < n do
          re(i) = buf(2 * i)
          im(i) = buf(2 * i + 1)
          i += 1
        true
      finally InterleavedScratch.release(n, buf)

  /** Power-of-two real FFT using JTransforms packed real layout. */
  def supportsReal(n: Int): Boolean =
    n >= 2 && (n & 1) == 0

  def realForward(re: Array[Double], im: Array[Double], onesided: Boolean): Boolean =
    val n = re.length
    if !supportsReal(n) || im.length != n then false
    else
      val half = n >> 1
      val buf = InterleavedScratch.borrow(half)
      try
        System.arraycopy(re, 0, buf, 0, n)
        plan(n).realForward(buf)
        re(0) = buf(0)
        im(0) = 0.0
        re(half) = buf(1)
        im(half) = 0.0
        var k = 1
        while k < half do
          val kr = buf(k << 1)
          val ki = buf((k << 1) + 1)
          re(k) = kr
          im(k) = ki
          if !onesided then
            re(n - k) = kr
            im(n - k) = -ki
          k += 1
        true
      finally InterleavedScratch.release(half, buf)

  /** Inverse of a onesided Hermitian spectrum, unscaled like [[pow2]]. */
  def realInverse(re: Array[Double], im: Array[Double]): Boolean =
    val n = re.length
    if !supportsReal(n) || im.length != n then false
    else
      val half = n >> 1
      val buf = InterleavedScratch.borrow(half)
      try
        buf(0) = re(0)
        buf(1) = re(half)
        var k = 1
        while k < half do
          buf(k << 1) = re(k)
          buf((k << 1) + 1) = im(k)
          k += 1
        plan(n).realInverse(buf, false)
        val unscaledFactor = if (n & (n - 1)) == 0 then 2.0 else 1.0
        var i = 0
        while i < n do
          re(i) = unscaledFactor * buf(i)
          im(i) = 0.0
          i += 1
        true
      finally InterleavedScratch.release(half, buf)

/** Length-`2n` interleaved complex scratch for JTransforms. */
private[internal] object InterleavedScratch:
  private val tl =
    new ThreadLocal[java.util.HashMap[Integer, java.util.ArrayDeque[Array[Double]]]]:
      override def initialValue() =
        new java.util.HashMap[Integer, java.util.ArrayDeque[Array[Double]]]()

  def borrow(n: Int): Array[Double] =
    val map = tl.get()
    val key = Integer.valueOf(n)
    var q = map.get(key)
    if q == null then
      q = new java.util.ArrayDeque[Array[Double]]()
      val _ = map.put(key, q)
    val hit = q.pollFirst()
    if hit != null then hit else new Array[Double](n * 2)

  def release(n: Int, buf: Array[Double]): Unit =
    val map = tl.get()
    val key = Integer.valueOf(n)
    var q = map.get(key)
    if q == null then
      q = new java.util.ArrayDeque[Array[Double]]()
      val _ = map.put(key, q)
    val _ = q.addFirst(buf)
