# signal4s 1.0 support matrix

Status relative to [`proposal.md`](../proposal.md) §6. This is an
implementation-snapshot matrix; each row has the boundary notes needed to
interpret “supported”.

## Required (1.0) — supported

| Area | Status | Module |
|---|---|---|
| Direct / FFT / OLA / circular convolution | supported; FFT/OLA require installed FFT hooks | `signal4s-core`, `signal4s-fft` |
| Correlation + lag axes | supported | `signal4s-core` |
| FIR / TF / SOS runners + state | supported | `signal4s-core` |
| Zero-phase (`filtfilt`-style) | supported for implemented padding paths; Gustafsson is not a complete implementation | `signal4s-core` |
| FFT / iFFT / rFFT | supported | `signal4s-fft` |
| Windows (periodic vs symmetric) | supported | `signal4s-fft` |
| Framing, STFT + inverse | supported for onesided analysis/synthesis; custom synthesis validation is incomplete | `signal4s-fft` |
| Periodogram / Welch | onesided supported; twosided estimation deferred | `signal4s-fft` |
| freqz / impulse / step / group delay | supported | `signal4s-design` |
| Mean / linear detrend | supported | `signal4s-fft` (Welch/periodogram) |
| Windowed-sinc FIR + Butterworth→SOS | supported | `signal4s-design` |
| `upfirdn` / rational polyphase resample | supported | `signal4s-core`, `signal4s-design` |
| Gale convolution operator + adjoint | supported for `Full`, `Valid`, and `Input(Zero)` | `signal4s-core` |
| Optional Ravel axis adapters | supported | `signal4s-ravel` |

## Explicitly deferred

Exotic windows; Chebyshev / elliptic / Bessel; Parks–McClellan; wavelets; AR/ME spectral estimators; peak finding; chirps / comms generators; audio codecs / I/O; image filters; Matlab / R / SciPy API compatibility layers; lazy signal graphs; Cats Effect / FS2 in core; autodiff; GPU; full state-space design conversions; band-stop Butterworth.

## Platforms

| Capability | JVM | Scala.js |
|---|---|---|
| Core numeric types + runners | yes | yes |
| FFT / STFT / Welch | yes | yes |
| Design (FIR / Butterworth) | yes | yes |
| SciPy fixture parity suites | yes | n/a (fixtures are JVM-checked) |
| Measured `Auto` cost model receipt | yes (portable) | model shared; platform-specific retune later |
| Ravel adapters | yes | yes |

## Benchmark receipts

Checked under [`benchmarks/receipts/`](../benchmarks/receipts/):

- `e3-direct-operator.md` — Gale convolution operator
- `e5-fft.md` — portable FFT
- `e6-auto-cost.md` — convolution `Auto` cost model
- `e6-fft-convolve.md` — FFT/OLA convolution
- `e10-sos.md` — SOS filter throughput
