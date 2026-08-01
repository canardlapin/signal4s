package signal4s.ravel

import gale.linalg.{DMat, Vec}
import signal4s.*
import signal4s.filter.Fir

class RavelSuite extends munit.FunSuite:

  test("AlongAxis maps FIR over rows with explicit copy"):
    val mat = DMat.tabulate(2, 4)((r, c) => (r + 1) * (c + 1).toDouble)
    val fir = Fir.causal(Vec(0.5, 0.5)).orThrow
    val out = AlongAxis
      .map(mat, Axis.Rows, ContiguousPolicy.CopyAlways)(x => fir.process(x))
      .orThrow
    assertEquals(out.rows, 2)
    assertEquals(out.cols, 4)
    // causal FIR [0.5,0.5] on [1,2,3,4] → [0.5, 1.5, 2.5, 3.5]
    assertEqualsDouble(out(0, 0), 0.5, 1e-12)
    assertEqualsDouble(out(0, 1), 1.5, 1e-12)
    assertEqualsDouble(out(1, 0), 1.0, 1e-12)

  test("ChannelBatch round-trip through matrix"):
    val batch = ChannelBatch
      .of(Vec(1.0, 2.0, 3.0), Vec(4.0, 5.0, 6.0))
      .orThrow
    val mat = ChannelBatch.toMatrix(batch, Axis.Rows).orThrow
    val back = ChannelBatch.fromMatrix(mat, Axis.Rows).orThrow
    assertEquals(back.channelCount, 2)
    assertEqualsDouble(back.channel(1)(2), 6.0, 0.0)

  test("RequireContiguous rejects strided column of row-major matrix"):
    val mat = DMat.tabulate(3, 3)((r, c) => (r * 3 + c).toDouble)
    // row-major contiguous rows OK
    val ok = AlongAxis.map(mat, Axis.Rows, ContiguousPolicy.RequireContiguous)(x => Right(x))
    assert(ok.isRight)
    // columns are strided in row-major layout
    val bad = AlongAxis.map(mat, Axis.Columns, ContiguousPolicy.RequireContiguous)(x => Right(x))
    assert(bad.isLeft)
