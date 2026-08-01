package signal4s.internal

import gale.linalg.DVec
import signal4s.Boundary

/** Boundary-aware sample reads for finite observations. */
private[signal4s] object SampleAccess:

  def apply(signal: DVec, index: Int, boundary: Boundary): Double =
    val n = signal.length
    if index >= 0 && index < n then signal(index)
    else
      boundary match
        case Boundary.Zero            => 0.0
        case Boundary.Constant(value) => value
        case Boundary.Clamp =>
          if n == 0 then 0.0
          else if index < 0 then signal(0)
          else signal(n - 1)
        case Boundary.Reflect =>
          if n == 0 then 0.0
          else if n == 1 then signal(0)
          else signal(reflectIndex(index, n))
        case Boundary.Symmetric =>
          if n == 0 then 0.0
          else signal(symmetricIndex(index, n))

  /** NumPy `reflect`: period 2*(n-1), edge not duplicated. */
  private def reflectIndex(index: Int, n: Int): Int =
    val period = 2 * (n - 1)
    var i = Math.floorMod(index, period)
    if i >= n then i = period - i
    i

  /** NumPy `symmetric`: period 2*n, edge duplicated. */
  private def symmetricIndex(index: Int, n: Int): Int =
    val period = 2 * n
    var i = Math.floorMod(index, period)
    if i >= n then i = period - 1 - i
    i
