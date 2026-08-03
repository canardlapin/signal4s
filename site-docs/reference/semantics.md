# Semantic conventions

The important interoperability decisions are explicit:

- `Signal` represents a finite observation on a regular lattice with a sample
  rate and start time.
- `Kernel` stores a zero-lag tap index. Causal and centered kernels therefore
  have different output coordinates even when their tap values are identical.
- `Full`, `Valid`, and `Input(boundary)` describe support and boundary meaning;
  `Direct`, `Fft`, and `OverlapAdd` describe computation method.
- `Input(Boundary.Zero)` is the signal-aligned equivalent of an odd centered
  “same” convolution. Other boundary policies have distinct meanings.
- Correlation returns an explicit lag axis; do not assume another library's
  lag sign convention.
- Window construction requires a `Periodic` or `Symmetric` convention.
- STFT frame alignment and phase shift are named choices, and real spectra are
  onesided with an explicit frequency axis.
- Welch and periodogram currently implement onesided output; twosided estimates
  are intentionally rejected at estimation time.
- Resampling is a stateful upsample–FIR–downsample operation. `flush` drains the
  FIR transient and does not add zero-insertion after the final input sample.

The repository's detailed numerical notes and fixture mappings remain in
[`docs/SEMANTICS.md`](https://github.com/canardlapin/signal4s/blob/main/docs/SEMANTICS.md)
and [`fixtures/CATALOG.md`](https://github.com/canardlapin/signal4s/blob/main/fixtures/CATALOG.md).

When porting a result from SciPy or gsignal, compare coordinates, support,
window convention, sides, and state/flush behavior before comparing numeric
values. See the [migration reference](migration.md).
