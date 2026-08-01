package signal4s.fft

import signal4s.ConvolutionMethod

/** Measured relative cost model for [[ConvolutionMethod.Auto]] (JVM portable).
  *
  * Coefficients come from `benchmarks/receipts/e6-auto-cost.md` (wall-clock
  * ratios of Direct / FFT / OLA on representative (N, M) pairs). Units are
  * arbitrary nanoseconds; only ordering matters.
  */
private[fft] object AutoCostModel:

  // Calibrated from benchmarks/receipts/e6-auto-cost.md (JVM portable).
  // On that matrix Direct won every cell; these keep modeled ordering aligned
  // until N·M grows large enough that FFT/OLA win on cost.
  private val cDirect: Double = 1.0
  private val cFft: Double = 12.0
  private val cOla: Double = 12.0

  /** Choose among Direct, FFT, and a default OLA block when beneficial. */
  def select(signalLength: Int, kernelLength: Int): ConvolutionMethod =
    val n = signalLength.toLong
    val m = kernelLength.toLong
    if n <= 0 || m <= 0 then ConvolutionMethod.Direct
    else
      val direct = cDirect * n * m
      val lFft = nextPow2(n + m - 1)
      val fft = cFft * lFft * log2(lFft)

      val block = defaultBlock(signalLength, kernelLength)
      val lOla = nextPow2(block.toLong + m - 1)
      val blocks = (n + block - 1) / block
      val ola = cOla * blocks * lOla * log2(lOla)

      if fft <= direct && fft <= ola then ConvolutionMethod.Fft
      else if ola < direct then ConvolutionMethod.OverlapAdd(block.toInt)
      else ConvolutionMethod.Direct

  /** Prefer blocks near a power-of-two FFT of modest size when N ≫ M. */
  def defaultBlock(signalLength: Int, kernelLength: Int): Int =
    val target = math.max(32, math.min(signalLength, 4 * kernelLength))
    math.max(1, nextPow2(target.toLong).toInt)

  private def nextPow2(min: Long): Long =
    var n = 1L
    while n < min do n <<= 1
    n

  private def log2(n: Long): Double =
    if n <= 1 then 1.0 else math.log(n.toDouble) / math.log(2.0)
