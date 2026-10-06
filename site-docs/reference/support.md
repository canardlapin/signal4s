# Support and maturity

This page describes the current implementation snapshot. It is not a promise
that every item has a stable published artifact or a hosted API reference.

## Implemented paths with boundaries

| Area | Current boundary |
|---|---|
| Convolution | Direct is available in core; FFT and overlap-add require the FFT hooks to be installed. |
| Correlation | Available with an explicit lag axis and signal4s lag convention. |
| FIR / transfer-function / SOS runners | Available with single-owner mutable state, snapshots, and reset. |
| Zero-phase filtering | Odd/constant/even padding paths exist; do not treat the named Gustafsson case as a complete implementation. |
| FFT / inverse FFT / real FFT | JVM and Scala.js implementations are tested; workspaces are single-owner. |
| STFT and inverse | Onesided analysis/synthesis is available; custom synthesis-window invertibility is not validated by construction. |
| Welch / periodogram | Onesided output is available; twosided estimation returns an explicit numerical failure. |
| FIR and Butterworth-to-SOS design | Available on JVM and Scala.js for the documented design families. |
| `upfirdn` / rational resampling | Batch and streaming paths are available; streaming uses owned checkpoints, restore/reset and explicit `flush`. Factors are GCD-reduced. |
| Gale convolution operator | Adjoint support is limited to `Full`, `Valid`, and `Input(Boundary.Zero)`. |
| Ravel axis adapters | Optional rows/columns/channel adapters are available on JVM and Scala.js. |

## Deferred or out of scope

Exotic windows; Chebyshev, elliptic, and Bessel designs; Parks–McClellan;
wavelets; AR/ME spectral estimators; peak finding; chirps and communications
generators; audio/image I/O; lazy graphs; Cats Effect/FS2 integration;
autodiff/GPU backends; full state-space conversions; and band-stop Butterworth.

## Platform and release evidence

The repository's JVM and Scala.js module suites are the stronger local
compatibility evidence. The fixture parity suite is JVM-only because it loads
committed SciPy fixtures. Performance receipts are machine- and JDK-specific;
read the [performance reference](performance.md) before using them as a
decision rule. The build currently reports `0.1.0-SNAPSHOT`, and no verified
Maven publication path is advertised here.

For a detailed matrix, see [`docs/SUPPORT.md`](https://github.com/canardlapin/signal4s/blob/main/docs/SUPPORT.md).
