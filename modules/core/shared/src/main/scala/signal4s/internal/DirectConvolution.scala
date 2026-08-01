package signal4s.internal

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*

/** Scalar convolution. Zero-boundary paths use scatter/gather kernels on raw
  * arrays; other boundaries use [[dotAt]] (the numeric oracle for laws).
  */
private[signal4s] object DirectConvolution:

  /** Full convolution array contents equal SciPy `convolve(..., mode="full")`
    * for the same tap sequence; [[Kernel.zeroLagIndex]] only shifts the time axis.
    *
    * Under [[Boundary.Zero]], `out[i] = Σ_j taps[j] · x[i - j]` independent of
    * the origin index (origin affects only the time axis metadata).
    */
  def full(signal: DVec, kernel: Kernel): DVec =
    fullZero(signal, kernel.taps)

  def valid(signal: DVec, kernel: Kernel): DVec =
    val n = signal.length
    val m = kernel.length
    val z = kernel.zeroLagIndex
    val first = m - 1 - z
    val last = n - 1 - z
    if last < first then DVec.zeros(0)
    else validZero(signal, kernel.taps, first, last - first + 1)

  def inputAligned(signal: DVec, kernel: Kernel, boundary: Boundary): DVec =
    boundary match
      case Boundary.Zero if kernel.zeroLagIndex == 0 =>
        inputZeroCausal(signal, kernel.taps)
      case Boundary.Zero =>
        inputZeroGeneral(signal, kernel)
      case other =>
        val n = signal.length
        val out = DVecBuilder.zeros(n)
        var i = 0
        while i < n do
          out(i) = dotAt(signal, kernel, i, other)
          i += 1
        out.result()

  private def fullZero(signal: DVec, taps: DVec): DVec =
    val n = signal.length
    val m = taps.length
    val outLen = n + m - 1
    val x = new Array[Double](n)
    val h = new Array[Double](m)
    signal.copyTo(x)
    taps.copyTo(h)
    val out = new Array[Double](outLen)
    PlatformDirect.scatterFull(x, h, out)
    adopt(out)

  private def validZero(signal: DVec, taps: DVec, first: Int, outLen: Int): DVec =
    val n = signal.length
    val m = taps.length
    if PlatformDirect.preferFullForValid(n, m) then
      val x = new Array[Double](n)
      val h = new Array[Double](m)
      signal.copyTo(x)
      taps.copyTo(h)
      val full = new Array[Double](n + m - 1)
      PlatformDirect.scatterFull(x, h, full)
      adoptRange(full, first, outLen)
    else
      val x = new Array[Double](n)
      val h = new Array[Double](m)
      signal.copyTo(x)
      taps.copyTo(h)
      val out = new Array[Double](outLen)
      var i = 0
      while i < outLen do
        val sampleIndex = first + i
        var acc = 0.0
        var j = 0
        // Fully interior: sampleIndex-j in [0,n) for all j in [0,m)
        val fullyInterior = sampleIndex >= m - 1 && sampleIndex < n
        if fullyInterior then
          while j < m do
            acc += h(j) * x(sampleIndex - j)
            j += 1
        else
          while j < m do
            val src = sampleIndex - j
            if src >= 0 && src < n then acc += h(j) * x(src)
            j += 1
        out(i) = acc
        i += 1
      adopt(out)

  /** Causal (`z = 0`): `y[i] = Σ_j h[j] x[i - j]`. */
  private def inputZeroCausal(signal: DVec, taps: DVec): DVec =
    val n = signal.length
    val m = taps.length
    if PlatformDirect.preferFullForCausalInput(n, m) then
      val x = new Array[Double](n)
      val h = new Array[Double](m)
      signal.copyTo(x)
      taps.copyTo(h)
      val full = new Array[Double](n + m - 1)
      PlatformDirect.scatterFull(x, h, full)
      adoptPrefix(full, n)
    else
      val x = new Array[Double](n)
      val h = new Array[Double](m)
      signal.copyTo(x)
      taps.copyTo(h)
      val out = new Array[Double](n)
      var i = 0
      while i < n do
        var acc = 0.0
        val jMax = math.min(m, i + 1)
        var j = 0
        while j < jMax do
          acc += h(j) * x(i - j)
          j += 1
        out(i) = acc
        i += 1
      adopt(out)

  /** General origin: `y[i] = Σ_j h[j] x[i - (j - z)]` with zero outside. */
  private def inputZeroGeneral(signal: DVec, kernel: Kernel): DVec =
    val n = signal.length
    val m = kernel.length
    val z = kernel.zeroLagIndex
    val x = new Array[Double](n)
    val h = new Array[Double](m)
    signal.copyTo(x)
    kernel.taps.copyTo(h)
    val out = new Array[Double](n)
    var i = 0
    while i < n do
      var acc = 0.0
      var j = 0
      while j < m do
        val src = i - j + z
        if src >= 0 && src < n then acc += h(j) * x(src)
        j += 1
      out(i) = acc
      i += 1
    adopt(out)

  private def adopt(values: Array[Double]): DVec =
    val builder = DVecBuilder.zeros(values.length)
    var i = 0
    while i < values.length do
      builder(i) = values(i)
      i += 1
    builder.result()

  private def adoptPrefix(values: Array[Double], length: Int): DVec =
    adoptRange(values, offset = 0, length)

  private def adoptRange(values: Array[Double], offset: Int, length: Int): DVec =
    val builder = DVecBuilder.zeros(length)
    var i = 0
    while i < length do
      builder(i) = values(offset + i)
      i += 1
    builder.result()

  def circular(signal: DVec, kernel: Kernel, period: Int): DVec =
    val out = DVecBuilder.zeros(period)
    val z = kernel.zeroLagIndex
    val m = kernel.length
    var i = 0
    while i < period do
      var acc = 0.0
      var j = 0
      while j < m do
        val lag = j - z
        val src = Math.floorMod(i - lag, period)
        val x =
          if src < signal.length then signal(src)
          else 0.0
        acc += kernel.taps(j) * x
        j += 1
      out(i) = acc
      i += 1
    out.result()

  /** y[n] = Σ_j taps[j] · x[n - (j - z)] */
  private def dotAt(
      signal: DVec,
      kernel: Kernel,
      n: Int,
      boundary: Boundary
  ): Double =
    val z = kernel.zeroLagIndex
    val m = kernel.length
    var acc = 0.0
    var j = 0
    while j < m do
      val lag = j - z
      acc += kernel.taps(j) * SampleAccess(signal, n - lag, boundary)
      j += 1
    acc

  def outputStart(signalStart: Seconds, kernel: Kernel, sampleRate: SampleRate): Seconds =
    Seconds.unsafe(signalStart.value + kernel.outputStartOffset(sampleRate).value)

  def validFirstSampleIndex(kernel: Kernel): Int =
    kernel.length - 1 - kernel.zeroLagIndex
