package signal4s.fft

import gale.linalg.DVec

/** Welch / periodogram estimate with retained diagnostics (no fabricated DoF). */
final case class WelchResult(
    frequencies: FrequencyAxis,
    power: DVec,
    scaling: SpectralScaling,
    sides: SpectralSides,
    segmentCount: Int,
    /** \(\sum_n w[n]^2\) used in SciPy-style normalization. */
    windowPower: Double,
    /** \(\left(\sum_n w[n]\right)^2\) used for spectrum scaling. */
    windowCoherentPower: Double,
    detrend: Detrend,
    average: AverageMethod,
    /** Optional equivalent degrees of freedom; absent unless a model is specified. */
    degreesOfFreedom: Option[Double],
    diagnostics: List[String]
)
