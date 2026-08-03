# Multichannel data with Ravel

The optional `signal4s-ravel` module adapts a one-dimensional operation to
matrix rows, columns, or planar channel batches. Core 1-D types do not depend
on Ravel.

```scala mdoc:silent
import gale.linalg.{DMat, Vec}
import signal4s.*
import signal4s.filter.Fir
import signal4s.ravel.*

val ravelFir = Fir.causal(Vec(0.25, 0.5, 0.25)).orThrow
val ravelMatrix = DMat.tabulate(4, 128)((channel, time) => math.sin(0.1 * time + channel))
val ravelFiltered =
  AlongAxis
    .map(ravelMatrix, Axis.Rows, ContiguousPolicy.CopyAlways)(row => ravelFir.process(row))
    .orThrow
assert(ravelFiltered.rows == ravelMatrix.rows)
assert(ravelFiltered.cols == ravelMatrix.cols)
```

`CopyAlways` gives the callback an owned contiguous vector. It is the portable
default for strided matrix views. `RequireContiguous` avoids the copy only when
the selected axis is already contiguous; otherwise it returns a
`NumericalFailure` rather than silently allocating.

The callback must return compatible lengths for every row or column. A length
change is reported as `SignalError.LengthMismatch`. If each channel needs
independent streaming history, create a runner per channel or use
`AlongAxis.mapChannels` with a stateful callback designed for that ownership
model; a single shared runner would mix channel state.

Next: [recover from typed errors](../reference/errors-and-recovery.md).
