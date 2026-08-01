package signal4s.laws

import gale.linalg.Vec
import signal4s.*
import signal4s.fft.*

class StftLawSuite extends munit.FunSuite:

  private val fs = SampleRate.hertz(8000.0).orThrow

  test("Hann 50% hop round-trip with edge frames"):
    val plan = StftPlan.hannPeriodic(64, fs, hop = Some(32)).orThrow
    StftLaws.dualComputedOnce(plan)
    val x = Vec.tabulate(200)(i => math.sin(0.1 * i) + 0.25 * math.cos(0.37 * i))
    StftLaws.roundTrip(plan, x)
    StftLaws.axesMatchPlan(plan, x)

  test("invalid hop rejected at plan construction"):
    val win = Window.fromSpec(WindowSpec.Hann(32, WindowConvention.Periodic)).orThrow
    assert(StftPlan(win, hop = 0, fs, nfft = 32).isLeft)
    assert(StftPlan(win, hop = 16, fs, nfft = 8).isLeft) // nfft < frame
