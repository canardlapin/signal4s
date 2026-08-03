# Core concepts

signal4s has three boundaries that are easy to blur in a conventional DSP
API. Keeping them separate makes numerical meaning and operational behavior
visible in types and values.

## 1. Coordinates are part of the value

`Signal` stores samples, a `SampleRate`, and the time of the first sample.
Operations that return a `Signal` calculate a new start coordinate; they do not
silently return a bare vector when alignment matters. Spectral results retain a
frequency axis, time-frequency results retain both axes, and correlations
retain a lag axis.

Some lower-level APIs intentionally return `DVec` because they operate on
unlabelled samples. When moving between those APIs and `Signal`, make the
sampling lattice explicit at the boundary.

## 2. A kernel has an origin

`Kernel.causal(taps)` places zero lag at tap `0`. `Kernel.at(taps, origin)` is
the general constructor. `Kernel.centeredOdd(taps)` places zero lag at the
middle tap and rejects even lengths, because an even-length “center” is a
choice that must be named rather than guessed.

The kernel origin determines the first coordinate of a `Full` convolution:
the output begins at `input.start - zeroLagIndex / sampleRate`.

## 3. Support and method are different

Select the returned support with `OutputRegion`; select the algorithm with
`ConvolutionMethod`. This lets a test compare direct and FFT execution without
changing a caller's requested coordinates, and lets `Auto` evolve without
changing a caller's semantic request.

The Gale convolution operator adds another contract: its currently implemented
adjoint is limited to `Full`, `Valid`, and `Input(Boundary.Zero)`. Ordinary
convolution and the operator are related, but they are not interchangeable
boundary APIs.

## 4. Descriptions, plans, runners, and workspaces

- A filter description (`Fir`, transfer function, or SOS cascade) is reusable
  and contains no streaming history.
- A runner owns mutable delay state. It is single-owner and not thread-safe;
  use `snapshot`, `restore`, or `reset` deliberately.
- An FFT or convolution plan is immutable and normally fixed to a length.
- An FFT workspace is mutable scratch owned by one caller. Allocate one per
  concurrent execution path rather than sharing it between threads.

Smart constructors return `Either[SignalError, A]`. Invalid rates, frequencies,
origins, lengths, unsupported regions, and unavailable numerical methods are
ordinary recoverable values until application code chooses to throw.

Next: [design a filter and inspect its spectrum](guides/filter-and-spectrum.md)
or [start with chunked execution](guides/streaming.md).
