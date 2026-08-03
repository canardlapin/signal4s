# Errors and recovery

Smart constructors and numerical operations return `Either[SignalError, A]`.
This makes invalid scientific configuration visible without forcing exceptions
on callers.

```scala mdoc:silent
import gale.linalg.Vec
import signal4s.*

val invalidKernel = Kernel.centeredOdd(Vec(1.0, 2.0, 3.0, 4.0))
assert(invalidKernel == Left(SignalError.KernelNotOddLength(4)))

val safeKernel =
  invalidKernel.fold(
    _ => Kernel.causal(Vec(1.0, 2.0, 1.0)).orThrow,
    identity
  )
assert(safeKernel.zeroLagIndex == 0)
```

Common recovery categories:

| Error | Typical response |
|---|---|
| `InvalidSampleRate`, `InvalidFrequency`, `FrequencyAboveNyquist` | Correct the physical configuration or report it to the caller. |
| `EmptySignal`, `EmptyKernel`, `InvalidInputLength` | Validate the batch or stream boundary before planning. |
| `KernelNotOddLength`, `InvalidKernelOrigin` | Choose and name the alignment explicitly with `Kernel.at`. |
| `LengthMismatch`, `InvalidFilterState` | Reuse a plan/state only with the length it was created for. |
| `UnsupportedOperatorRegion` | Use zero-extension operator regions until adjoints for other boundaries exist. |
| `NumericalFailure` | Inspect the operation/detail pair; this includes unavailable FFT methods and deferred spectral or edge methods. |

Use pattern matching when the caller can recover differently by category. Use
`.orThrow` only at an application boundary where failure is intentionally fatal,
or in small examples where it keeps the semantic point visible.

Next: [read the semantic conventions](semantics.md) before comparing results
with another DSP library.
