package signal4s.multirate

import gale.linalg.{DVec, Vec}
import signal4s.*

class MultirateSuite extends munit.FunSuite:

  private val h = Vec(0.25, 0.5, 0.25)
  private val x = Vec.tabulate(8)(i => (i + 1).toDouble)

  test("RateRatio GCD reduction"):
    val r = RateRatio(4, 6).orThrow
    assertEquals(r.up, 2)
    assertEquals(r.down, 3)
    assert(RateRatio(0, 1).isLeft)
    assert(RateRatio(1, -1).isLeft)

  test("upfirdn up=1 down=1 matches full FIR convolution length"):
    val y = Upfirdn(h, x, 1, 1).orThrow
    assertEquals(y.length, x.length + h.length - 1)

  test("upfirdn validates empty kernels and gives empty input a defined result"):
    assert(Upfirdn(DVec.zeros(0), x, 1, 1).isLeft)
    assertClose(Upfirdn(h, DVec.zeros(0), 3, 2).orThrow, DVec.zeros(0))
    assertEquals(Upfirdn.outputLength(h.length, 0, 3, 2), 0)

  test("upfirdn known SciPy vectors"):
    assertClose(Upfirdn(h, x, 1, 1).orThrow, Vec(0.25, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 5.75, 2.0))
    assertClose(Upfirdn(h, x, 1, 2).orThrow, Vec(0.25, 2.0, 4.0, 6.0, 5.75))
    assertClose(Upfirdn(h, x, 2, 3).orThrow, Vec(0.25, 1.0, 1.75, 2.5, 3.25, 4.0))

  test("general polyphase recurrence matches materialized upsample FIR"):
    val signal = Vec.tabulate(19)(i => math.sin(0.17 * i) + 0.05 * (i % 4))
    val taps = Vec.tabulate(17)(i => math.cos(0.23 * i) / (i + 1.0))
    val ratios =
      for
        up <- 2 to 7
        down <- 2 to 7
        if gcd(up, down) == 1
      yield (up, down)

    ratios.foreach { case (up, down) =>
      val actual = Upfirdn(taps, signal, up, down).orThrow
      val full = Array.ofDim[Double]((signal.length - 1) * up + taps.length)
      var xi = 0
      while xi < signal.length do
        var k = 0
        while k < taps.length do
          full(xi * up + k) += signal(xi) * taps(k)
          k += 1
        xi += 1
      val expected = DVec.tabulate(Upfirdn.outputLength(taps.length, signal.length, up, down)) { oi =>
        full(oi * down)
      }
      assertClose(actual, expected, tol = 1e-12)
    }

  test("streaming consume+flush ≡ batch upfirdn"):
    List((1, 1), (2, 1), (1, 2), (2, 3), (3, 2)).foreach { case (up, down) =>
      val batch = Upfirdn(h, x, up, down).orThrow
      val r = PolyphaseResampler(h, up, down).orThrow
      val a = r.consume(x.slice(0, 3).copy).orThrow
      val b = r.consume(x.slice(3, 8).copy).orThrow
      val c = r.flush().orThrow
      val stream = concat(a, b, c)
      assertClose(stream, batch, tol = 1e-12)
    }

  test("identity ratio with unit impulse FIR is identity"):
    val delta = Vec(1.0)
    val y = Upfirdn(delta, x, RateRatio.identity).orThrow
    assertClose(y, x)

  test("Resample convenience operations are thin, equivalent upfirdn contracts"):
    val ratio = RateRatio(6, 4).orThrow
    assertEquals((ratio.up, ratio.down), (3, 2))
    assertEquals(ratio.invert, RateRatio(2, 3).orThrow)
    assertEqualsDouble(ratio.scale, 1.5, 1e-15)
    assert(!ratio.isIdentity)
    assert(RateRatio.identity.isIdentity)

    assertClose(Resample(x, h, ratio).orThrow, Upfirdn(h, x, ratio).orThrow)
    assertClose(Resample(x, h, 3, 2).orThrow, Upfirdn(h, x, 3, 2).orThrow)
    assertClose(Resample.decimate(x, h, 2).orThrow, Upfirdn(h, x, 1, 2).orThrow)
    assertClose(Resample.interpolate(x, h, 3).orThrow, Upfirdn(h, x, 3, 1).orThrow)
    assertClose(Resample.identityFiltered(x, h).orThrow, Upfirdn(h, x, 1, 1).orThrow)
    assert(Resample.decimate(x, h, 0).isLeft)
    assert(Resample.interpolate(x, h, -1).isLeft)

  test("polyphase bank lengths cover prototype"):
    val bank = PolyphaseBank.decompose(h, 2).orThrow
    assertEquals(bank.nPhases, 2)
    assertEquals(bank.phases(0).length + bank.phases(1).length, h.length)

  test("polyphase banks retain every residue class and validate degenerate inputs"):
    val prototype = Vec(1.0, 2.0, 3.0, 4.0, 5.0)
    val bank = PolyphaseBank.decompose(prototype, 3).orThrow
    assertEquals(bank.prototypeLength, prototype.length)
    assertEquals(bank.maxPhaseLength, 2)
    assertClose(bank.phases(0), Vec(1.0, 4.0))
    assertClose(bank.phases(1), Vec(2.0, 5.0))
    assertClose(bank.phases(2), Vec(3.0))
    val sparse = PolyphaseBank.decompose(prototype, 7).orThrow
    assertEquals(sparse.nPhases, 7)
    assertEquals(sparse.maxPhaseLength, 1)
    assertEquals(sparse.phases(6).length, 0)
    assert(PolyphaseBank.decompose(DVec.zeros(0), 2).isLeft)
    assert(PolyphaseBank.decompose(prototype, 0).isLeft)

  test("streaming resampler flush is idempotent and closes the input"):
    assert(PolyphaseResampler(DVec.zeros(0), RateRatio.identity).isLeft)
    val resampler = PolyphaseResampler(h, up = 2, down = 3).orThrow
    assertEquals(resampler.samplesConsumed, 0)
    assertClose(resampler.consume(DVec.zeros(0)).orThrow, DVec.zeros(0))
    assertClose(resampler.flush().orThrow, DVec.zeros(0))
    assert(resampler.isFlushed)
    assertClose(resampler.flush().orThrow, DVec.zeros(0))
    assert(resampler.consume(Vec(1.0)).isLeft)

  private def concat(parts: DVec*): DVec =
    val n = parts.map(_.length).sum
    val out = gale.linalg.DVecBuilder.zeros(n)
    var o = 0
    parts.foreach { p =>
      var i = 0
      while i < p.length do
        out(o) = p(i)
        o += 1
        i += 1
    }
    out.result()

  private def gcd(a: Int, b: Int): Int =
    var x = a
    var y = b
    while y != 0 do
      val r = x % y
      x = y
      y = r
    x

  private def assertClose(a: gale.linalg.DVec, b: gale.linalg.DVec, tol: Double = 1e-12): Unit =
    assertEquals(a.length, b.length)
    var i = 0
    while i < a.length do
      assert(math.abs(a(i) - b(i)) <= tol, s"i=$i ${a(i)} vs ${b(i)}")
      i += 1
