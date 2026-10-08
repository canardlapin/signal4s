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
Its `samplesConsumed` counter is `Long`, so long streams preserve global phase
past the Int boundary. This is a source API widening from the earlier Int counter.
The runner clocks inserted zeros and owns one prototype-sized delay line;
output vectors are owned allocations. `snapshot` copies an immutable in-memory
checkpoint; `restore` accepts only the same reduced ratio and exact prototype.
Snapshots retain global clocks, registers and flush state. `reset` explicitly
starts a fresh stream. Serialized checkpoint formats and acquisition-segment
identity remain consumer concerns.
`flush` drains the FIR transient and marks the resampler flushed; consuming
after that returns an error, while a second flush returns an empty vector.
There is no implicit filter flush for `FirRunner` or `SosRunner`: they process
exactly the samples supplied.

Do not share a runner or workspace concurrently. If two execution paths need
the same filter description, give each path its own runner.

Next: [reuse immutable FFT plans](fft-and-plans.md) or [map a filter across channels](multichannel.md).

## FIR finite refusal and capacities

Causal FIR construction requires finite coefficients. A runner stages delay state
and output before committing a block: nonfinite input/state or overflowing output/
delay arithmetic returns an error and preserves the completed prefix and caller
output destination. Restored/initial state must also be finite. Ordinary successful
DF-II transposed arithmetic and partition conventions are unchanged. Rounded
underflow follows binary64; exact cancellation of overflowing intermediate terms
is not claimed.

```scala mdoc:silent
val finiteFirGuide = signal4s.filter.Fir.causal(gale.linalg.Vec(0.25, 0.5, 0.25)).orThrow
val finiteFirResources = finiteFirGuide.resources
assert(finiteFirResources.stateBytes == 16)
assert(finiteFirResources.coefficientBytes == 32)
assert(finiteFirResources.processIntoScratchBytes(5).orThrow == 56)
assert(finiteFirResources.ownedProcessAdditionalBytes(5).orThrow == 136)
```

Resources describe primitive array payloads: taps plus the shared unit denominator,
delay registers, independent snapshots, staged block/state and owned-output copies.
Caller input/destination, stack scalars, object/reference/allocator/GC/RSS overhead
and escaped output lifetime are excluded. Inspection allocates no runner/workspace.
These are capacity facts, not a zero-allocation or complete-workload performance
claim. Runners/destinations remain single-owner, not thread-safe.
