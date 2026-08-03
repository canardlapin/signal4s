# Migrating from SciPy or gsignal

`signal4s` is not an API port. Map meaning first, then choose a method.

| SciPy / gsignal idea | signal4s starting point |
|---|---|
| `convolve(..., mode="full")` | `Convolution(x, Kernel.causal(h), OutputRegion.Full)` |
| `convolve(..., mode="valid")` | `OutputRegion.Valid` |
| odd-kernel `mode="same"` | `Kernel.centeredOdd(h)` with `OutputRegion.Input(Boundary.Zero)` |
| `lfilter(b, a, x)` | a `Fir`, `DigitalTransferFunction`, or `SecondOrderCascade` runner |
| `sosfilt(sos, x)` | `SecondOrderCascade.fromSosMatrix(...).process(x)` |
| `filtfilt(..., padtype="odd")` | `ZeroPhase.filter(..., EdgeTreatment.OddPad)` |
| `get_window(..., fftbins=...)` | `WindowSpec.*` with `Periodic` or `Symmetric` |
| `welch` / `periodogram` | `WelchPlan` / `Periodogram` with onesided output |
| `firwin` | `FirDesign.lowPass`, `highPass`, or `bandPass` |
| `butter(..., output="sos")` | `Butterworth.lowPass(...).sos` or another implemented design |
| `upfirdn` / `resample_poly` | `Upfirdn` / `ResamplePoly` |
| along-axis multi-channel processing | `signal4s.ravel.AlongAxis` |

The [filter and spectrum guide](../guides/filter-and-spectrum.md) shows a
complete design-to-Welch path. The [semantic reference](semantics.md) covers
the cases where matching names do not imply matching coordinates or edge
behavior.

The detailed repository migration table and examples remain in
[`docs/MIGRATION.md`](https://github.com/canardlapin/signal4s/blob/main/docs/MIGRATION.md).
