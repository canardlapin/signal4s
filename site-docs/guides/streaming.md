# Streaming, state, and flush

Batch convenience methods create a fresh runner. For chunked input, create one
runner and retain ownership of it across calls.

```scala mdoc:silent
import gale.linalg.Vec
import signal4s.*
import signal4s.filter.*
import signal4s.multirate.*

val streamingFir = Fir.causal(Vec(0.25, 0.5, 0.25)).orThrow
val streamingRunner = streamingFir.newRunner()
val streamingFirst = streamingRunner.process(Vec(1.0, 2.0, 3.0)).orThrow
val streamingCheckpoint = streamingRunner.snapshot
val streamingSecond = streamingRunner.process(Vec(4.0, 5.0)).orThrow

streamingRunner.restore(streamingCheckpoint).orThrow
val streamingReplay = streamingRunner.process(Vec(4.0, 5.0)).orThrow
assert(streamingSecond.length == streamingReplay.length)
assert(
  (0 until streamingSecond.length).forall { i =>
    math.abs(streamingSecond(i) - streamingReplay(i)) < 1e-12
  }
)

val streamingResampler =
  PolyphaseResampler(Vec(1.0, 0.5, 0.25), up = 2, down = 3).orThrow
val streamingPart = streamingResampler.consume(Vec(1.0, 2.0, 3.0)).orThrow
val streamingTail = streamingResampler.flush().orThrow
assert(streamingPart.length + streamingTail.length > 0)
```

FIR, transfer-function, and SOS runners are single-owner mutable state. Use
`processInto` when the caller owns an output buffer, `reset` to return to zero
state, and `snapshot`/`restore` for checkpointing or deterministic replay.

`PolyphaseResampler.consume` may be called repeatedly until end of stream.
`flush` drains the FIR transient and marks the resampler flushed; consuming
after that returns an error, while a second flush returns an empty vector.
There is no implicit filter flush for `FirRunner` or `SosRunner`: they process
exactly the samples supplied.

Do not share a runner or workspace concurrently. If two execution paths need
the same filter description, give each path its own runner.

Next: [reuse immutable FFT plans](fft-and-plans.md) or [map a filter across channels](multichannel.md).
