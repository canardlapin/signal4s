# Design a filter and estimate a spectrum

This is the common analysis path for a sampled measurement: design a stable
Butterworth filter as second-order sections, process the samples, then estimate
a Welch power spectrum.

```scala mdoc:silent
import gale.linalg.Vec
import signal4s.*
import signal4s.design.*
import signal4s.fft.*

val filterGuideRate = SampleRate.hertz(1000.0).orThrow
val filterGuideInput =
  Vec.tabulate(256)(i => math.sin(2 * math.Pi * 40 * i / filterGuideRate.hertz))

val filterGuideIir =
  Butterworth.lowPass(4, Frequency.hertz(100.0).orThrow, filterGuideRate).orThrow
val filterGuideOutput = filterGuideIir.sos.process(filterGuideInput).orThrow

val filterGuideWindow =
  Window.fromSpec(WindowSpec.Hann(64, WindowConvention.Periodic)).orThrow
val filterGuideWelch =
  WelchPlan(
    filterGuideWindow,
    hop = 32,
    sampleRate = filterGuideRate,
    nfft = 64,
    detrend = Detrend.Mean
  ).orThrow
val filterGuideSpectrum = filterGuideWelch.estimate(filterGuideOutput).orThrow

assert(filterGuideSpectrum.frequencies.length == 33)
assert(filterGuideSpectrum.power.length == 33)
```

The SOS representation avoids the numerical sensitivity of expanding a high-
order IIR into one polynomial transfer function. The window convention is
explicit: `Periodic` is the usual FFT-bin convention, while `Symmetric` is a
different finite-window definition.

The current Welch implementation is onesided. A `Twosided` plan can be
constructed but estimation returns a `NumericalFailure`; choose the default
`SpectralSides.Onesided` until that capability is implemented. The same
limitation applies to `Periodogram`.

For zero-phase filtering, use an implemented edge treatment such as
`EdgeTreatment.OddPad`. The enum also names SciPy's `Gustafsson` method, but
the current implementation does not provide a true Gustafsson calculation for
all filter forms; see [support and maturity](../reference/support.md).

Next: [carry a filter across chunks](streaming.md) or [reuse FFT plans](fft-and-plans.md).

## Owned bounded periodograms

`BoundedPeriodogramPlan` estimates exactly one full frame using a caller-owned
workspace. Its portable radix-2/Bluestein route owns FFT buffers and tables;
it does not depend on native scratch pools or global numerical caches. Odd,
prime and padded FFT lengths are admitted up to the explicit shape cap.
Inspect capacities before reserving a workspace. Primitive payload capacities
include both direction tables, chirps/kernel arrays, frame scratch and owned
output. Object/reference/allocator/GC/RSS costs are excluded; this is not a
heap limit or a complete-workload benchmark.

```scala mdoc:silent
val boundedGuidePlan = BoundedPeriodogramPlan(
  filterGuideWindow, filterGuideRate, 65, Detrend.Mean,
  SpectralScaling.Density
).orThrow
val boundedGuideCapacity = boundedGuidePlan.resources
assert(boundedGuideCapacity.workspaceBytes > 0)
val boundedGuideWorkspace = boundedGuidePlan.newWorkspace()
val boundedGuideFrame = filterGuideOutput.slice(0, 64)
val boundedGuideFirst = boundedGuidePlan.estimateInto(boundedGuideFrame, boundedGuideWorkspace).orThrow
boundedGuidePlan.estimateInto(Vec.zeros(64), boundedGuideWorkspace).orThrow
assert(boundedGuideFirst.power.toSeq.exists(_ > 0))
assert(boundedGuideFirst.segmentCount == 1)
assert(boundedGuideFirst.degreesOfFreedom.isEmpty)
```

A workspace belongs to its exact plan and is single-owner/not thread-safe.
Returned power vectors own their storage. Foreign workspaces, shapes, nonfinite
inputs or unrepresentable power return errors; no partial power result escapes.
Mean detrending uses anchored differences and exact normalized sums, preserving
small variation under a large constant offset. Frame magnitudes are normalized
before FFT and restored with exponent arithmetic after spectral normalization.
Welch's mean accumulation uses exact sums/normalized readout so an overflowing
sum of finite frame powers does not invalidate a finite mean. Rounded underflow
follows binary64. Anchored differences/detrended magnitudes exceeding finite
capacity are explicit failures; this is not exact real-arithmetic detrending.
The legacy raw-median convention remains unchanged, without a SciPy median-bias
correction; only its numerical midpoint and size admission are improved.

These APIs estimate linear power. Logarithmic baseline correction and its
averaging order belong to a separately declared operation.


## Bounded centered resampling

`ResamplePolyPlan` designs one immutable Kaiser5 prototype for the reduced rational
rate. `inputWindow` computes the exact observed support of an absolute output
window; `processWindow` consumes precisely that support and retains segment-global
phase. No source-length array or replay from ordinal0 is required. Outside the
segment is zero padding; symmetric prototype delay is removed, as in SciPy
resample_poly with constant padding. Output length is ceil(inputLength*up/down).
Long clocks and support use exact integer arithmetic, with typed capacity refusal.
Finite input and late overflow expose no reusable partial result. ResamplePoly
batch convenience shares this numeric core. Raw Upfirdn/streaming resampler remain
separately described APIs, without this plan's bounded/finite contract.

```scala mdoc
import signal4s.design.ResamplePolyPlan
val centeredPlan = ResamplePolyPlan(2, 3).orThrow
val requiredSupport = centeredPlan.inputWindow(1000L, 301L, 8).orThrow
val observedSupport = gale.linalg.Vec.tabulate(requiredSupport.count)(_ => 1.0)
val centeredWindow = centeredPlan.processWindow(observedSupport, requiredSupport.start, 1000L, 301L, 8).orThrow
```

`coefficientBytes`, `maximumInputCount`, `inputHalo` and
`processAdditionalBytes` describe retained/staged primitive numeric payloads.
The output estimate conservatively includes builder and completed owned storage.
Caller input, metadata/objects/BigInt/allocator/GC/RSS and plan-design construction
are excluded. These are engineering capacities, not end-to-end performance claims.


Kaiser taps use the positive I0 series through beta500 and a scaled/logarithmic
large-argument expansion beyond it, avoiding overflowing normalization. The formulas
follow [NIST DLMF10.25.2](https://dlmf.nist.gov/10.25.E2) and
[NIST DLMF10.40.1](https://dlmf.nist.gov/10.40.E1). Independent frozen SciPy i0e ratios
cover symmetric/periodic windows across beta0..1000; extreme finite beta also stays
finite. Underflowed tails may round to zero. Matching resample_poly smoke tolerance
is1e-12 relative plus1e-12 absolute. No new special-function library is required.
