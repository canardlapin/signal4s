# signal4s

`signal4s` is a Scala 3 library for finite, regularly sampled one-dimensional
signals. It combines convolution, filtering, Fourier analysis, spectral
estimation, and rational resampling while keeping coordinates and algorithm
choices explicit.

This documentation is generated from the `site-docs/` source tree. Its Scala
examples are compiled by `sbt docs/tlSite` against the JVM modules. The shared
library is cross-built for the JVM and Scala.js; the site gate is a JVM guide
gate, not a Scala.js proof.

## Learning path

1. [Getting started](getting-started.md) — construct a signal, convolve it,
   and inspect the output coordinates.
2. [Concepts](concepts.md) — understand axes, kernel origin, support, method,
   and ownership.
3. [Filter and spectrum guide](guides/filter-and-spectrum.md) — design a
   Butterworth SOS filter and estimate a Welch spectrum.
4. [Streaming guide](guides/streaming.md) — process chunks, checkpoint state,
   and flush a rational resampler.
5. [FFT and plans](guides/fft-and-plans.md) — install the FFT backend and reuse
   immutable plans with single-owner workspaces.
6. [Multichannel guide](guides/multichannel.md) — apply 1-D operations along a
   matrix axis with the optional Ravel module.

## Reference

- [Errors and recovery](reference/errors-and-recovery.md)
- [Semantic conventions](reference/semantics.md)
- [Support and maturity](reference/support.md)
- [Migration from SciPy / gsignal](reference/migration.md)
- [Performance interpretation](reference/performance.md)

Generated Scaladoc is not published yet. For a local API reference, run
`sbt coreJVM/doc`, `sbt fftJVM/doc`, or the corresponding module task.

## Scope

The current checkout is a `0.1.0-SNAPSHOT` implementation snapshot. The
support page distinguishes implemented paths from deliberately deferred
features and known limitations; the site should not be read as a release or
artifact-availability claim.
