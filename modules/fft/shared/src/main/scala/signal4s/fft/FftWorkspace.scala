package signal4s.fft

/** Single-owner mutable scratch for an [[FftPlan]] / [[RealFftPlan]].
  *
  * Workspaces must not be cached inside plans. Not thread-safe.
  */
final class FftWorkspace private[fft] (
    private[fft] val re: Array[Double],
    private[fft] val im: Array[Double],
    private[fft] val scratchRe: Array[Double],
    private[fft] val scratchIm: Array[Double]
):
  def length: Int = re.length

  private[fft] def clear(): Unit =
    var i = 0
    while i < re.length do
      re(i) = 0.0
      im(i) = 0.0
      i += 1

object FftWorkspace:
  private[fft] def apply(
      re: Array[Double],
      im: Array[Double],
      scratchRe: Array[Double],
      scratchIm: Array[Double]
  ): FftWorkspace =
    new FftWorkspace(re, im, scratchRe, scratchIm)
