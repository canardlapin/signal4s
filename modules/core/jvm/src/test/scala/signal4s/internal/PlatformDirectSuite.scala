package signal4s.internal

/** Platform acceleration must retain scalar full-convolution semantics. */
class PlatformDirectSuite extends munit.FunSuite:

  test("macOS selects Accelerate for dense supported FIRs"):
    val isMac = sys.props.getOrElse("os.name", "").toLowerCase.contains("mac")
    if isMac then
      assertEquals(
        PlatformDirect.backendFor(n = 257, m = 17),
        "accelerate",
        clue = PlatformDirect.accelerateDiagnostic
      )
    else
      assertNotEquals(PlatformDirect.backendFor(n = 257, m = 17), "accelerate")

  test("accelerated-or-scalar scatter matches reference on a dense FIR"):
    assertScatterMatches(n = 257, m = 17)

  test("dispatch thresholds preserve scalar FIR semantics"):
    assertEquals(PlatformDirect.backendFor(n = 17, m = 8), "short")
    assertEquals(PlatformDirect.backendFor(n = 17, m = 9), "unrolled")
    assertEquals(PlatformDirect.backendFor(n = 64, m = 9), "blocked")
    assertNotEquals(PlatformDirect.backendFor(n = 64, m = 2045), "accelerate")

    List((17, 8), (17, 9), (64, 9), (64, 16), (64, 2044), (64, 2045)).foreach {
      case (n, m) => assertScatterMatches(n, m)
    }

  test("causal input borrows full scatter only from native or vector tiers"):
    val backend = PlatformDirect.backendFor(n = 1024, m = 48)
    assertEquals(
      PlatformDirect.preferFullForCausalInput(n = 1024, m = 48),
      backend == "accelerate" || backend == "vector"
    )
    assertEquals(
      PlatformDirect.preferFullForValid(n = 1024, m = 48),
      PlatformDirect.preferFullForCausalInput(n = 1024, m = 48)
    )
    assert(!PlatformDirect.preferFullForCausalInput(n = 48, m = 1024))

  private def assertScatterMatches(n: Int, m: Int): Unit =
    val x = Array.tabulate(n) { i =>
      if i % 7 == 0 then 0.0 else math.sin(0.07 * i) + 0.01 * (i % 5)
    }
    val h = Array.tabulate(m) { i =>
      if i % 11 == 0 then 0.0 else math.cos(0.13 * i) / (i + 1.0)
    }
    val got = Array.ofDim[Double](n + m - 1)
    val expected = Array.ofDim[Double](n + m - 1)

    PlatformDirect.scatterFull(x, h, got)
    var i = 0
    while i < n do
      var j = 0
      while j < m do
        expected(i + j) += x(i) * h(j)
        j += 1
      i += 1

    i = 0
    while i < got.length do
      assertEqualsDouble(got(i), expected(i), 1e-11, clue = s"n=$n m=$m i=$i")
      i += 1
