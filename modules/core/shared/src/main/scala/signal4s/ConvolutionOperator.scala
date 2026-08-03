package signal4s

import gale.linalg.{DVec, DoubleLinearOperator, LinAlgError, MutableDVec}
import signal4s.internal.DirectConvolution

/** Finite batch convolution as a Gale `DoubleLinearOperator`.
  *
  * Only zero-extension regions are supported: [[OutputRegion.Full]],
  * [[OutputRegion.Valid]], and [[OutputRegion.Input]] with [[Boundary.Zero]].
  * Reflective / clamped / constant / symmetric boundaries are rejected because
  * their adjoints accumulate multiple extended positions onto one sample and are
  * not yet implemented.
  *
  * This is a finite batch map. Streaming filter runners must not extend
  * `DoubleLinearOperator`.
  */
final class ConvolutionOperator private (
    val kernel: Kernel,
    val inputLength: Int,
    val region: OutputRegion
) extends DoubleLinearOperator:

  override val cols: Int = inputLength
  override val rows: Int = ConvolutionOperator.outputLength(inputLength, kernel, region)

  override def applyTo(x: DVec, into: MutableDVec): Unit =
    if x.length != cols then throw LinAlgError.VectorLengthMismatch(cols, x.length)
    if into.length != rows then throw LinAlgError.VectorLengthMismatch(rows, into.length)
    val y = region match
      case OutputRegion.Full              => DirectConvolution.full(x, kernel)
      case OutputRegion.Valid             => DirectConvolution.valid(x, kernel)
      case OutputRegion.Input(Boundary.Zero) =>
        DirectConvolution.inputAligned(x, kernel, Boundary.Zero)
      case other =>
        throw LinAlgError.UnsupportedOperation(
          s"convolution operator forward unsupported for $other"
        )
    into := y

  override def transposeApplyTo(y: DVec, into: MutableDVec): Unit =
    if y.length != rows then throw LinAlgError.VectorLengthMismatch(rows, y.length)
    if into.length != cols then throw LinAlgError.VectorLengthMismatch(cols, into.length)
    into.clear()
    region match
      case OutputRegion.Full =>
        adjointFull(y, into)
      case OutputRegion.Valid =>
        adjointValid(y, into)
      case OutputRegion.Input(Boundary.Zero) =>
        adjointInputZero(y, into)
      case other =>
        throw LinAlgError.UnsupportedOperation(
          s"convolution operator adjoint unsupported for $other"
        )

  /** (A* y)[k] = Σ_j h[j] y[k+j] for Full / SciPy-style indexing. */
  private def adjointFull(y: DVec, into: MutableDVec): Unit =
    val m = kernel.length
    val outLen = y.length
    var k = 0
    while k < cols do
      var acc = 0.0
      var j = 0
      while j < m do
        val yi = k + j
        if yi >= 0 && yi < outLen then acc += kernel.taps(j) * y(yi)
        j += 1
      into(k) = acc
      k += 1

  /** Adjoint of Valid under zero extension. */
  private def adjointValid(y: DVec, into: MutableDVec): Unit =
    val m = kernel.length
    val first = DirectConvolution.validFirstSampleIndex(kernel)
    var k = 0
    while k < cols do
      var acc = 0.0
      var j = 0
      while j < m do
        // `first` already contains the kernel origin. Subtracting the origin
        // again shifts the adjoint of centered valid convolution.
        val yi = k - first + j
        if yi >= 0 && yi < y.length then acc += kernel.taps(j) * y(yi)
        j += 1
      into(k) = acc
      k += 1

  /** (A* y)[k] = Σ_j h[j] y[k+j-z] for Input(Zero). */
  private def adjointInputZero(y: DVec, into: MutableDVec): Unit =
    val m = kernel.length
    val z = kernel.zeroLagIndex
    var k = 0
    while k < cols do
      var acc = 0.0
      var j = 0
      while j < m do
        val yi = k + j - z
        if yi >= 0 && yi < y.length then acc += kernel.taps(j) * y(yi)
        j += 1
      into(k) = acc
      k += 1

object ConvolutionOperator:

  private[signal4s] def make(
      kernel: Kernel,
      inputLength: Int,
      region: OutputRegion
  ): ConvolutionOperator =
    new ConvolutionOperator(kernel, inputLength, region)

  private[signal4s] def outputLength(
      inputLength: Int,
      kernel: Kernel,
      region: OutputRegion
  ): Int =
    region match
      case OutputRegion.Full  => inputLength + kernel.length - 1
      case OutputRegion.Valid =>
        math.max(0, inputLength - kernel.length + 1)
      case OutputRegion.Input(_) => inputLength
