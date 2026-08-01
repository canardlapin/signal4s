package signal4s.laws

import gale.linalg.Vec
import signal4s.*

class KernelLawSuite extends munit.ScalaCheckSuite:
  test("causal kernel lag range"):
    val kernel = Kernel.causal(Vec(1.0, 2.0, 3.0)).orThrow
    KernelLaws.lagRangeAgreesWithOrigin(kernel)
    KernelLaws.causalOriginIsZero(kernel)

  test("centeredOdd kernel lag range"):
    val kernel = Kernel.centeredOdd(Vec(0.25, 0.5, 0.25)).orThrow
    KernelLaws.lagRangeAgreesWithOrigin(kernel)
    KernelLaws.centeredOddOriginIsMiddle(kernel)

  test("explicit origin lag range"):
    val kernel = Kernel.at(Vec(1.0, 2.0, 3.0, 4.0), zeroLagIndex = 2).orThrow
    KernelLaws.lagRangeAgreesWithOrigin(kernel)
    assertEquals(kernel.minLag, -2)
    assertEquals(kernel.maxLag, 1)
