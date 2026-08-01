# E3 direct convolution operator bench

generated: 2026-08-01T03:15:28.645975Z
inputLength: 10000
kernelLength: 64
region: Full
iters: 100
totalMs: 67.249
msPerCall: 0.6725
approxBytesPerCall: 242787
sink: 0.0

Notes: allocating `op(x)` path (immutable result). Destination-writing
forms arrive with FFT/OLA plans; this receipt baselines the current
DirectConvolution-backed operator.
