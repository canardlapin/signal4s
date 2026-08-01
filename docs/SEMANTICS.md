# Intentional semantic differences

These are deliberate. Do not “fix” them to match SciPy silently.

## Kernel origin

- SciPy `convolve` treats tap arrays as causal sequences; “same” centering is a separate mode.
- signal4s `Kernel` always carries `zeroLagIndex`. Even-length kernels cannot use `centeredOdd`.
- Full convolution **sample values** match SciPy for the same tap sequence; origin shifts the **time axis**, not the array contents of `OutputRegion.Full`.

## `OutputRegion.Input` vs SciPy `same`

- Named `Input`, not `Same`.
- SciPy `mode="same"` with an odd kernel maps to `Kernel.centeredOdd` + `Input(Boundary.Zero)`.
- Other boundaries (`Reflect`, `Symmetric`, …) are first-class and not SciPy mode strings.

## Correlation lag sign

\[
r_{xy}[\ell] = \sum_n x[n]\, y[n+\ell]
\]

SciPy `correlate` uses the opposite lag convention for the same formula presentation. Fixtures document the signal4s convention explicitly (`smoke.correlate_raw_impulse_odd`).

## Frequencies are typed

- Cutoffs are `Frequency` + `SampleRate`, never bare `Double`.
- Values above Nyquist are rejected (`FrequencyAboveNyquist`).

## Window convention is mandatory

- `WindowConvention.Periodic` ↔ SciPy `fftbins=True`
- `WindowConvention.Symmetric` ↔ SciPy `fftbins=False`
- No default that silently picks one.

## STFT phase

- signal4s matches SciPy `ShortTimeFFT(..., phase_shift=None)` (raw windowed rFFT).
- SciPy’s default `phase_shift=0` applies a centering phase we do not apply.

## Welch degrees of freedom

- `WelchResult.degreesOfFreedom` is `None` unless a model is supplied.
- We do not fabricate DoF from segment count alone.

## SOS pairing

- Butterworth **ZPK** matches SciPy.
- SOS **section order** from `Sos.fromZerosPolesGain` may differ from `zpk2sos`; frequency / impulse responses match.

## Streaming resampler flush

- `PolyphaseResampler.consume` + `flush` matches batch `upfirdn`.
- Zeros are inserted **between** input samples only, not after the last sample, before the FIR drain (SciPy `upfirdn` xu length).

## Representations stay distinct

- TF, ZPK, and SOS are different values. Prefer SOS for high-order IIR execution.
- Design reports carry warnings (pole magnitude, high order) instead of expanding `SignalError`.
