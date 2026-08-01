package signal4s.laws

import gale.linalg.Vec
import signal4s.*
import signal4s.multirate.*

class MultirateLawSuite extends munit.FunSuite:

  test("ratio reduction preserves meaning"):
    MultirateLaws.ratioReductionPreservesMeaning(12, 18)
    MultirateLaws.ratioReductionPreservesMeaning(7, 5)

  test("identity with delta FIR"):
    MultirateLaws.identityDeltaIsIdentity(Vec.tabulate(16)(_.toDouble))

  test("stream matches batch for several ratios"):
    val h = Vec(0.2, 0.3, 0.3, 0.2)
    val x = Vec.tabulate(20)(i => math.cos(0.3 * i))
    List((1, 1), (2, 1), (1, 2), (3, 2)).foreach { case (u, d) =>
      MultirateLaws.streamMatchesBatch(h, x, u, d)
    }

  test("chunk latency bounded before flush"):
    val h = Vec(0.25, 0.5, 0.25)
    MultirateLaws.chunkLatencyBounded(h, RateRatio(2, 3).orThrow, chunk = 8)
