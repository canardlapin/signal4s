package signal4s.internal

/** Scala.js direct FIR kernel must match the scalar definition, including sparse inputs. */
class PlatformDirectSuite extends munit.FunSuite:

  test("scatterFull matches scalar reference across unrolled and tail shapes"):
    List((17, 7), (64, 17), (65, 18)).foreach { case (n, m) =>
      val x = Array.tabulate(n) { i =>
        if i % 6 == 0 then 0.0 else math.sin(0.19 * i) + 0.02 * (i % 3)
      }
      val h = Array.tabulate(m) { i =>
        if i % 5 == 0 then 0.0 else math.cos(0.23 * i) / (i + 1.0)
      }
      val actual = Array.ofDim[Double](n + m - 1)
      val expected = Array.ofDim[Double](n + m - 1)

      PlatformDirect.scatterFull(x, h, actual)
      var i = 0
      while i < n do
        var j = 0
        while j < m do
          expected(i + j) += x(i) * h(j)
          j += 1
        i += 1

      i = 0
      while i < actual.length do
        assertEqualsDouble(actual(i), expected(i), 1e-12, clue = s"n=$n m=$m i=$i")
        i += 1
    }
