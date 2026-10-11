# Owned one-frame complex STFT analysis

`BoundedStftPlan` analyzes one complete real observation frame and returns
detached positive-frequency complex coefficients. The immutable plan retains
original and unit-normalized window taps, sample rate, FFT size, detrend and
normalization. It has no hop, recording length, clock-origin assumption, inverse
transform or overlap-add condition. Scientific callers place frames on their
own global grid and carry validity/units/timeline authority separately.

`inspect(frameLength,nfft)` reports primitive/limb payload capacities without
allocating frame buffers, FFT arrays or tables. `withinBudget` refuses insufficient
construction/run capacity before copying window taps. `newWorkspace(budget)`
admits the retained plan/frame/result peak before allocating single-owner FFT
scratch. Keep one workspace per concurrent worker. Objects/references, allocator,
GC/RSS, caller input buffers and escaped retained results have separate budgets.

The result's phase origin is the extracted frame's first sample. The coefficient
convention is the unshifted rFFT (`ShortTimeFFT phase_shift=None`); there is no
implicit centering phase or interior-bin doubling. Explicit modes are:

- Raw: windowed Fourier sum, without an amplitude divisor.
- Magnitude: divide by abs(sum(window)), preserving the window's signed phase.
- Density: divide by sqrt(sampleRate*sum(window²)).

These are complex normalizations. Folded one-sided power requires a separately
declared interior-bin factor; complex coefficients themselves are not doubled.
The caller must distinguish input units, input units/√Hz, squared magnitude and
folded density. Frequency bins use the retained sample rate and nfft.

None/Mean/Linear detrending precedes windowing/FFT padding and shares owning
SpectralFrame arithmetic with periodograms. Only the transform is padded after
a complete frame; no incomplete observations are silently filled or normalized.
Wrong input length, nonfinite input, foreign workspace, undefined normalization,
unrepresentable frequency spacing or normalized coefficient refuse. A raw zero
window is valid analysis; magnitude/density require their nonzero gains.

Unit-normalized window taps and exponent restoration avoid raw amplitude/power
overflow and tiny-window gain underflow before normalization. Representability
is still binary64, not exact real arithmetic; anchored detrend residual limits
remain explicit and underflow rounds according to binary64.

The portable owned radix2/Bluestein engine and scratch layout are shared with
BoundedPeriodogramPlan. There are no native/global FFT scratch pools. Detached
columns stay unchanged after workspace reuse; failures produce no reusable
column. Frame streaming, scientific masks/segments, persistence and cancellation
are consumer capabilities, not inferred from this numerical primitive.
