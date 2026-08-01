package signal4s.internal

/** JVM hot kernels for zero-boundary direct convolution. */
private[internal] object PlatformDirect:

  private val Block = 64

  private val isArm: Boolean =
    val arch = sys.props.getOrElse("os.arch", "").toLowerCase
    arch.contains("aarch64") || arch.contains("arm64") || arch == "arm"

  private val isMac: Boolean =
    sys.props.getOrElse("os.name", "").toLowerCase.contains("mac")

  /** One whole-FIR vDSP call. It is macOS-only and retains the portable paths
    * whenever Accelerate/native access is unavailable.
    */
  private val accelerateScatter: Option[(Array[Double], Array[Double], Array[Double]) => Unit] =
    if !isMac then None
    else
      try
        val cls = Class.forName("signal4s.internal.AccelerateDirect")
        val available = cls.getDeclaredMethod("available").invoke(null).asInstanceOf[Boolean]
        if available then
          val m = cls.getDeclaredMethod(
            "convolveFull",
            classOf[Array[Double]],
            classOf[Array[Double]],
            classOf[Array[Double]]
          )
          Some((x, h, out) => { val _ = m.invoke(null, x, h, out); () })
        else None
      catch case _: Throwable => None

  /** Vector API FIR when hardware lanes ≥ 4 (typically x86 AVX). Skipped on
    * aarch64 where preferred width is 2 and wider species are emulated.
    */
  private val vectorScatter: Option[(Array[Double], Array[Double], Array[Double]) => Unit] =
    if isArm then None
    else
      try
        val cls = Class.forName("signal4s.internal.VectorScatter")
        val worthwhile = cls.getDeclaredMethod("worthwhile").invoke(null).asInstanceOf[Boolean]
        if worthwhile then
          val m = cls.getDeclaredMethod(
            "scatterFull",
            classOf[Array[Double]],
            classOf[Array[Double]],
            classOf[Array[Double]]
          )
          Some((x, h, out) => { val _ = m.invoke(null, x, h, out); () })
        else None
      catch case _: Throwable => None

  /** `out[i+j] += x[i] * h[j]` for all i,j (full FIR scatter). */
  def scatterFull(x: Array[Double], h: Array[Double], out: Array[Double]): Unit =
    val n = x.length
    val m = h.length
    backendFor(n, m) match
      case "accelerate" =>
        // The predicate behind `backendFor` proves the option is nonempty.
        accelerateScatter.get(x, h, out)
      case "vector" =>
        vectorScatter.get(x, h, out)
      case "short" =>
        scatterFullShort(x, h, out, n, m)
      case "blocked" =>
        scatterFullBlocked(x, h, out, n, m)
      case _ =>
        scatterFullUnrolled(x, h, out, n, m)

  /** Selected kernel for a shape; package-visible only for backend parity tests. */
  private[internal] def backendFor(n: Int, m: Int): String =
    if m <= 8 then "short"
    else if accelerateScatter.nonEmpty && m >= 16 && m <= 2044 && n >= 64 then "accelerate"
    else if vectorScatter.nonEmpty && m >= 16 && n >= 64 then "vector"
    else if n >= Block then "blocked"
    else "unrolled"

  /** Full scatter is faster than a scalar causal gather only on native/vector tiers.
    *
    * The full tail is intentionally avoided when the kernel exceeds the signal:
    * causal output then needs materially less work than a complete convolution.
    */
  private[internal] def preferFullForCausalInput(n: Int, m: Int): Boolean =
    m <= n &&
      (backendFor(n, m) match
        case "accelerate" | "vector" => true
        case _                        => false)

  private[internal] def preferFullForValid(n: Int, m: Int): Boolean =
    preferFullForCausalInput(n, m)

  /** Native-backend diagnostic for macOS test failures and benchmark receipts. */
  private[internal] def accelerateDiagnostic: String =
    if !isMac then "not macOS"
    else
      try
        Class
          .forName("signal4s.internal.AccelerateDirect")
          .getDeclaredMethod("diagnostic")
          .invoke(null)
          .asInstanceOf[String]
      catch case error: Throwable => s"diagnostic unavailable: ${error.getClass.getSimpleName}"

  private def scatterFullShort(
      x: Array[Double],
      h: Array[Double],
      out: Array[Double],
      n: Int,
      m: Int
  ): Unit =
    var i = 0
    while i < n do
      val xi = x(i)
      if xi != 0.0 then
        var j = 0
        while j < m do
          out(i + j) += xi * h(j)
          j += 1
      i += 1

  private def scatterFullUnrolled(
      x: Array[Double],
      h: Array[Double],
      out: Array[Double],
      n: Int,
      m: Int
  ): Unit =
    val m8 = m & ~7
    var i = 0
    while i < n do
      val xi = x(i)
      if xi != 0.0 then
        var j = 0
        while j < m8 do
          val base = i + j
          out(base) += xi * h(j)
          out(base + 1) += xi * h(j + 1)
          out(base + 2) += xi * h(j + 2)
          out(base + 3) += xi * h(j + 3)
          out(base + 4) += xi * h(j + 4)
          out(base + 5) += xi * h(j + 5)
          out(base + 6) += xi * h(j + 6)
          out(base + 7) += xi * h(j + 7)
          j += 8
        while j < m do
          out(i + j) += xi * h(j)
          j += 1
      i += 1

  private def scatterFullBlocked(
      x: Array[Double],
      h: Array[Double],
      out: Array[Double],
      n: Int,
      m: Int
  ): Unit =
    var i0 = 0
    while i0 < n do
      val i1 = math.min(i0 + Block, n)
      var j = 0
      while j < m do
        val hj = h(j)
        if hj != 0.0 then
          var i = i0
          val i4 = i1 - ((i1 - i0) & 3)
          while i < i4 do
            out(i + j) += hj * x(i)
            out(i + j + 1) += hj * x(i + 1)
            out(i + j + 2) += hj * x(i + 2)
            out(i + j + 3) += hj * x(i + 3)
            i += 4
          while i < i1 do
            out(i + j) += hj * x(i)
            i += 1
        j += 1
      i0 += Block
