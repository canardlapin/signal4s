package signal4s.laws

import gale.linalg.{DVec, DVecBuilder}
import signal4s.*
import signal4s.filter.*
import munit.Assertions

object FilterLaws extends Assertions:

  def assertClose(actual: DVec, expected: DVec, tol: Double = 1e-9): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      val scale = math.max(1.0, math.max(math.abs(actual(i)), math.abs(expected(i))))
      assert(
        math.abs(actual(i) - expected(i)) <= tol * scale,
        s"mismatch at $i: ${actual(i)} vs ${expected(i)}"
      )
      i += 1

  /** FIR from-rest batch equals causal Input(Zero) convolution. */
  def firBatchEqualsConvolution(fir: Fir, x: DVec): Unit =
    val yFir = fir.process(x).orThrow
    val yConv =
      Convolution(x, fir.kernel, OutputRegion.Input(Boundary.Zero)).orThrow
    assertClose(yFir, yConv, tol = 1e-12)

  /** run(x ++ y, s0) ≈ run(y, s_x) after running x from s0. */
  def chunkConcatenation(
      processChunk: (DVec, FilterState) => (DVec, FilterState),
      zero: FilterState,
      x: DVec,
      y: DVec
  ): Unit =
    val xy = concat(x, y)
    val (full, _) = processChunk(xy, zero)
    val (yx, sx) = processChunk(x, zero)
    val (yy, _) = processChunk(y, sx)
    val chunked = concat(yx, yy)
    assertClose(full, chunked, tol = 1e-9)

  def snapshotRestart(
      process: (DVec, FilterState) => (DVec, FilterState),
      zero: FilterState,
      x: DVec,
      y: DVec
  ): Unit =
    val (_, mid) = process(x, zero)
    val (cont, _) = process(y, mid)
    val (_, mid2) = process(x, zero)
    assertEquals(mid.length, mid2.length)
    var i = 0
    while i < mid.length do
      assertEqualsDouble(mid.values(i), mid2.values(i), 1e-15)
      i += 1
    val (cont2, _) = process(y, mid2)
    assertClose(cont, cont2, tol = 1e-12)

  def independentRunners(make: () => FirRunner, x: DVec): Unit =
    val r1 = make()
    val r2 = make()
    val y1 = r1.process(x).orThrow
    val y2a = r2.process(x).orThrow
    assertClose(y1, y2a)
    // Advance r2 only; reset r1 and confirm it still matches a fresh pass.
    val _ = r2.process(x).orThrow
    r1.reset()
    val y1b = r1.process(x).orThrow
    assertClose(y1b, y2a)

  private def concat(a: DVec, b: DVec): DVec =
    val out = DVecBuilder.zeros(a.length + b.length)
    var i = 0
    while i < a.length do
      out(i) = a(i)
      i += 1
    var j = 0
    while j < b.length do
      out(a.length + j) = b(j)
      j += 1
    out.result()
