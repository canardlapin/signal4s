package signal4s.fft

import signal4s.*

class WindowSuite extends munit.FunSuite:

  test("Hann periodic gains are finite and positive"):
    val w = Window.fromSpec(WindowSpec.Hann(32, WindowConvention.Periodic)).orThrow
    assert(w.coherentGain > 0.0)
    assert(w.powerGain > 0.0)
    assert(w.enbwBins > 0.0)
    assertEquals(w.length, 32)

  test("convention is required and retained"):
    val s = Window.fromSpec(WindowSpec.Hamming(8, WindowConvention.Symmetric)).orThrow
    val p = Window.fromSpec(WindowSpec.Hamming(8, WindowConvention.Periodic)).orThrow
    assertEquals(s.convention, WindowConvention.Symmetric)
    assertEquals(p.convention, WindowConvention.Periodic)
    assert(math.abs(s.taps(1) - p.taps(1)) > 1e-6)
