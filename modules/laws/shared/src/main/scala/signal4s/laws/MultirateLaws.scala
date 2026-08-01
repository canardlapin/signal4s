package signal4s.laws

import gale.linalg.DVec
import signal4s.*
import signal4s.multirate.*
import munit.Assertions

object MultirateLaws extends Assertions:

  def ratioReductionPreservesMeaning(up: Int, down: Int): Unit =
    val r = RateRatio(up, down).orThrow
    assertEquals(r.up.toLong * down.toLong, r.down.toLong * up.toLong)
    assertEquals(RateRatio(r.up, r.down).orThrow, r)

  def identityDeltaIsIdentity(signal: DVec): Unit =
    val y = Upfirdn(DVec.fromSeq(Seq(1.0)), signal, RateRatio.identity).orThrow
    assertEquals(y.length, signal.length)
    var i = 0
    while i < signal.length do
      assertEqualsDouble(y(i), signal(i), 0.0)
      i += 1

  def streamMatchesBatch(h: DVec, signal: DVec, up: Int, down: Int): Unit =
    val batch = Upfirdn(h, signal, up, down).orThrow
    val r = PolyphaseResampler(h, up, down).orThrow
    val mid = math.max(1, signal.length / 2)
    val a = r.consume(signal.slice(0, mid).copy).orThrow
    val b = r.consume(signal.slice(mid, signal.length).copy).orThrow
    val c = r.flush().orThrow
    val n = a.length + b.length + c.length
    assertEquals(n, batch.length)
    var i = 0
    def check(part: DVec): Unit =
      var j = 0
      while j < part.length do
        assertEqualsDouble(part(j), batch(i), 1e-12)
        i += 1
        j += 1
    check(a); check(b); check(c)

  def chunkLatencyBounded(h: DVec, ratio: RateRatio, chunk: Int): Unit =
    val r = PolyphaseResampler(h, ratio).orThrow
    val x = DVec.tabulate(chunk)(i => math.sin(0.2 * i))
    val y = r.consume(x).orThrow
    // Before flush, output cannot exceed the batch length for this prefix
    val maxBatch = Upfirdn.outputLength(h.length, chunk, ratio.up, ratio.down)
    assert(y.length <= maxBatch)
