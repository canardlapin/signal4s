package signal4s

import signal4s.testkit.SciPyFixture

/** JVM-only smoke check that committed SciPy fixtures are present and well-formed.
  * Does not invoke Python.
  */
class FixtureSmokeSuite extends munit.FunSuite:

  private val requiredIds = List(
    "smoke.impulse",
    "smoke.step",
    "smoke.sine_on_bin",
    "smoke.sine_between_bins",
    "smoke.convolve_full_causal_odd",
    "smoke.convolve_full_even_kernel_origin0",
    "smoke.convolve_valid_causal_odd",
    "smoke.convolve_same_odd_as_centered_input_zero",
    "smoke.correlate_raw_impulse_odd",
    "smoke.fft_forward_sine_on_bin"
  )

  test("all committed smoke fixtures parse with pinned metadata and tolerances"):
    val fixtures = SciPyFixture.allSmoke()
    assert(fixtures.nonEmpty, "no smoke fixtures found")
    val ids = fixtures.map(_.id).toSet
    assert(requiredIds.forall(ids.contains), clue = s"missing ids: ${requiredIds.filterNot(ids.contains)}")
