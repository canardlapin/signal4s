# Bounded cross spectra and Morlet coefficients

These FFT-module estimators accept an explicit primitive-payload budget before
allocating their workspaces. The budget excludes caller inputs, object headers,
garbage collection and process RSS. Returned arrays own their storage.

## Common-window cross-spectral density

`BoundedCsdPlan` estimates `mean(conj(X_i) * X_j)` with one-sided density
scaling. Every channel contributes to the same complete frames. Supply explicit
window taps, FFT length, sample rate and detrending; the caller is responsible
for excluding frames that cross acquisition gaps or contain missing samples.
A refused frame leaves the accumulator's accepted count and sums unchanged.

```scala mdoc:silent
import gale.linalg.Vec
import signal4s.*
import signal4s.fft.*

val spectralBudget = BigInt(64) * 1024 * 1024
val spectralRate = SampleRate.hertz(128).orThrow
val spectralWindow = Window.fromSpec(
  WindowSpec.Rectangular(16, WindowConvention.Periodic)).orThrow
val crossPlan = BoundedCsdPlan.withinBudget(
  2, spectralWindow, spectralRate, 16, spectralBudget).orThrow
val crossState = crossPlan.newAccumulator(spectralBudget).orThrow
val tone = Vec.tabulate(16)(i => math.cos(2 * math.Pi * i / 16))
crossState.accumulateFrame(Vector(tone, tone)).orThrow
val crossResult = crossState.finishOwned().orThrow.get
assert(crossResult.acceptedWindows == 1)
assert(crossResult(0, 1, 1).imaginary == 0.0)
```

The result stores the full Hermitian matrix in
`((row * channels + column) * bins + frequency)` order. Mixed channel units
produce their product per Hz. DC and even-length Nyquist bins are not doubled.
Exact sums average the finite binary64 products; they do not make FFTs or
products exact. Empty accumulators return `None`. Merge only independent
accumulators constructed from the same plan instance. Each accumulator has one
owner; simultaneous mutations are unsupported.

`crossResult.averageBins(groups, budget)` optionally computes arithmetic means
over nonempty, ordered, disjoint groups of original FFT-bin indices. Groups may
omit bins; every matrix entry uses the same groups. `sourceFrequencies` and
`frequencyGroups` retain the original axis and exact bin membership. The result
remains a density, not integrated band power. Admission includes the retained
source matrix and the new output builders; inspect it with
`CsdGroupedMatrix.inspect(channels, bins, groups)` before allocating the result.

The portable reference tests compare four pinned SciPy cases. The reviewed
error bound is `1e-30 + 1e-10 * sqrt(p_i * p_j)`, with each `p` the maximum
reference diagonal over all bins. Original per-bin bounds failed at spectral
nulls; their diagnostics are retained. Hermitian symmetry and a real diagonal
are exact construction properties. No coherence, independence, inferential
degrees of freedom or performance claim is supplied.

## Finite-support Morlet coefficients

`BoundedMorletPlan` follows the MNE convention: a Gaussian complex oscillation
truncated at five temporal standard deviations and normalized to discrete L2
norm `sqrt(2)`. The optional analytic zero-mean correction is explicit.
Coefficients keep the input unit; squared magnitudes are power in its square,
not a spectral density or calibrated analytic-signal amplitude.

```scala mdoc:silent
val waveletPlan = BoundedMorletPlan.withinBudget(
  spectralRate, Vector(8.0, 16.0, 32.0), Vector(3.0, 5.0, 7.0),
  spectralBudget, zeroMean = true).orThrow
val waveletInput = Vec.tabulate(256)(i => if i == 128 then 1.0 else 0.0)
val waveletBlock = waveletPlan.analyzeBlock(
  waveletInput, firstCenter = 0, samples = 256, spectralBudget).orThrow
assert(!waveletBlock.isValid(0, 0))
assert(waveletBlock.isValid(0, 128))
```

Storage is frequency-by-time, with separate real, imaginary and validity
components. A coefficient requires its entire frequency-specific support.
Invalid edges contain zero placeholders with no scientific value. For tiled
input, include the required halos and request only the tile's center samples.
The caller must separately apply acquisition validity and requested-run
boundaries; this API receives one finite contiguous vector, not missing-data
metadata. It never bridges separate vectors or pads their edges.

The portable tests compare twelve pinned MNE kernel/convolution cases, both
zero-mean policies, haloed chunks, ownership and refusal behavior. Direct
convolution is bounded by the supplied output budget and a maximum of 65,535
kernel samples. No decimation, inverse transform or speed comparison is offered.
