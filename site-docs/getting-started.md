# Getting started

## Build from a checkout

This snapshot does not have verified published dependency coordinates. Clone
the repository and run the build from its root. Gale is resolved from an
explicit `-Dsignal4s.gale.build=...` override, a sibling `../gale` checkout,
or the pinned revision in `build.sbt`.

```bash
sbt coreJVM/test fftJVM/test designJVM/test ravelJVM/test lawsJVM/test
sbt coreJS/test fftJS/test designJS/test ravelJS/test lawsJS/test
```

## Preserve coordinates through convolution

`Signal` couples samples to a regular sampling lattice. A `Kernel` carries an
explicit zero-lag tap, so the output origin is a consequence of the operation
rather than an undocumented convention.

```scala mdoc:silent
import gale.linalg.Vec
import signal4s.*

val gettingStartedRate = SampleRate.hertz(1000.0).orThrow
val gettingStartedInput =
  Signal(Vec(1.0, 2.0, 3.0, 4.0), Sampling(gettingStartedRate)).orThrow
val gettingStartedKernel = Kernel.centeredOdd(Vec(0.25, 0.5, 0.25)).orThrow
val gettingStartedOutput =
  Convolution(gettingStartedInput, gettingStartedKernel, OutputRegion.Full).orThrow

assert(gettingStartedOutput.length == 6)
assert(gettingStartedOutput.start == Seconds.of(-0.001).orThrow)
```

The three output-region meanings are distinct:

| Region | Coordinates and support |
|---|---|
| `Full` | Every finite-support overlap under zero extension; the coordinate range can extend before and after the input. |
| `Valid` | Only positions where every tap overlaps an observed input sample. |
| `Input(boundary)` | Exactly the input coordinates, with an explicit boundary policy. |

`ConvolutionMethod` is a separate choice. `Direct`, `Fft`, and `OverlapAdd`
change how the result is computed, not what `Full`, `Valid`, or `Input` means.
Use `.orThrow` only in short examples; application code should normally match
on the returned `Either[SignalError, A]`.

Next: read [the core concepts](concepts.md), then choose a task-oriented guide.
