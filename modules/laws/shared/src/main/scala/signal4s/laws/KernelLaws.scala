package signal4s.laws

import signal4s.*
import munit.Assertions

/** Reusable kernel laws expressed against the public API. */
object KernelLaws extends Assertions:

  /** Lag range agrees with origin and tap count. */
  def lagRangeAgreesWithOrigin(kernel: Kernel): Unit =
    assertEquals(kernel.minLag, -kernel.zeroLagIndex)
    assertEquals(kernel.maxLag, kernel.length - 1 - kernel.zeroLagIndex)
    assertEquals(kernel.maxLag - kernel.minLag + 1, kernel.length)

  /** Causal kernels place lag zero at the first tap. */
  def causalOriginIsZero(kernel: Kernel): Unit =
    assertEquals(kernel.zeroLagIndex, 0)
    assertEquals(kernel.minLag, 0)

  /** Odd-centered kernels place lag zero at the middle tap. */
  def centeredOddOriginIsMiddle(kernel: Kernel): Unit =
    assertEquals(kernel.length % 2, 1, "centeredOdd kernels must be odd length")
    assertEquals(kernel.zeroLagIndex, kernel.length / 2)
    assertEquals(kernel.minLag, -kernel.maxLag)
