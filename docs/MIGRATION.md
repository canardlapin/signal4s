# Migrating from SciPy / gsignal

`signal4s` is **not** an API port. Map meaning first, then pick a method.

## Quick map

| SciPy / gsignal | signal4s |
|---|---|
| `convolve(x, h, mode='full')` | `Convolution(x, Kernel.causal(h), OutputRegion.Full)` |
| `convolve(..., mode='valid')` | `OutputRegion.Valid` |
| `convolve(..., mode='same')` (odd `h`) | `Kernel.centeredOdd(h)` + `OutputRegion.Input(Boundary.Zero)` |
| `correlate` (SciPy lag sign) | `Correlate` — see [SEMANTICS.md](SEMANTICS.md); lag convention differs |
| `lfilter(b, a, x)` | `Fir` / `DigitalTransferFunction` / `SecondOrderCascade` runners |
| `sosfilt(sos, x)` | `SecondOrderCascade.fromSosMatrix(...).process(x)` |
| `filtfilt(..., padtype='odd')` | `ZeroPhase.filter(..., EdgeTreatment.OddPad)` |
| `get_window(w, n, fftbins=)` | `WindowSpec.*(n, Periodic\|Symmetric)` — convention **required** |
| `ShortTimeFFT` / `stft` | `StftPlan` — use `phase_shift=None` mapping (raw windowed rFFT) |
| `welch` / `periodogram` | `WelchPlan` / `Periodogram` |
| `firwin` | `FirDesign.lowPass` / `highPass` / `bandPass` |
| `butter` → SOS | `Butterworth.lowPass` → `DesignedIir.sos` |
| `upfirdn` | `signal4s.multirate.Upfirdn` |
| `resample_poly` | `signal4s.design.ResamplePoly` |
| Multi-channel along axis | `signal4s.ravel.AlongAxis` (optional module) |

## Example: SciPy-style smooth + spectrum

```scala
import gale.linalg.Vec
import signal4s.*
import signal4s.design.*
import signal4s.fft.*

val fs = SampleRate.hertz(1000.0).orThrow
val x = Vec.tabulate(256)(i => math.sin(2 * math.Pi * 40 * i / fs.hertz))

// Design + execute SOS (prefer SOS over expanded TF)
val iir = Butterworth.lowPass(4, Frequency.hertz(100.0).orThrow, fs).orThrow
val y = iir.sos.process(x).orThrow

// Welch density
val win = Window.fromSpec(WindowSpec.Hann(64, WindowConvention.Periodic)).orThrow
val welch = WelchPlan(win, hop = 32, sampleRate = fs, nfft = 64, detrend = Detrend.Mean).orThrow
val psd = welch.estimate(y).orThrow
val _ = psd.power
```

## Example: gsignal-style FIR convolve

gsignal often centers kernels implicitly. In signal4s the origin is always explicit:

```scala
import gale.linalg.Vec
import signal4s.*

val h = Kernel.centeredOdd(Vec(0.25, 0.5, 0.25)).orThrow
val x = Vec(1.0, 2.0, 3.0, 4.0)
val same = Convolution(x, h, OutputRegion.Input(Boundary.Zero)).orThrow
val _ = same
```

## Multi-channel (Ravel)

```scala
import gale.linalg.{DMat, Vec}
import signal4s.*
import signal4s.filter.Fir
import signal4s.ravel.*

val fir = Fir.causal(Vec(0.25, 0.5, 0.25)).orThrow
val mat = DMat.tabulate(4, 128)((ch, t) => math.sin(0.1 * t + ch))
val filtered =
  AlongAxis.map(mat, Axis.Rows, ContiguousPolicy.CopyAlways)(row => fir.process(row)).orThrow
val _ = filtered
```

Core 1-D types never import `signal4s-ravel`.
