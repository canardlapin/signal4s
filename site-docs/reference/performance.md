# Performance interpretation

Performance evidence in this repository is scoped evidence, not a universal
ranking. The committed receipts record a particular machine, JDK, platform,
input grid, warmup, and trial protocol. They should answer “what was measured?”
before they answer “which method is fastest?”

The current comparison material includes direct versus FFT/overlap-add
convolution, the `Auto` cost model, SOS throughput, and `upfirdn`. The E11
performance work is still active, so receipt ranges may change and should not
be copied into a release announcement without rerunning the benchmark gate.

Use `ConvolutionMethod.Direct` when predictable small-input behavior matters;
use `Fft` or `OverlapAdd` when the input/kernel regime and installed backend
justify it; use `Auto` when the measured cost model is acceptable for the target
platform. Reuse plans and give concurrent paths separate workspaces.

Detailed notes live in [`docs/PERFORMANCE.md`](https://github.com/canardlapin/signal4s/blob/main/docs/PERFORMANCE.md),
and raw receipts live under [`benchmarks/receipts/`](https://github.com/canardlapin/signal4s/tree/main/benchmarks/receipts).
