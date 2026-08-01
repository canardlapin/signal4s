#!/usr/bin/env python3
"""Generate committed smoke fixtures. Not invoked by CI."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scipy
from scipy import signal

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "smoke"


def versions() -> dict:
    return {
        "scipy": scipy.__version__,
        "numpy": np.__version__,
        "gsignal": None,
    }


def write(fixture: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{fixture['id']}.json"
    path.write_text(json.dumps(fixture, indent=2, sort_keys=True) + "\n")
    print(f"wrote {path.relative_to(ROOT)}")


def arr(a: np.ndarray) -> list[float]:
    return np.asarray(a, dtype=np.float64).tolist()


def main() -> None:
    n = 32
    impulse = np.zeros(n)
    impulse[0] = 1.0
    step = np.ones(n)
    fs = 1000.0
    t = np.arange(n) / fs
    sine_on_bin = np.sin(2.0 * np.pi * (4.0 * fs / n) * t)
    sine_between = np.sin(2.0 * np.pi * 12.3 * t)

    write(
        {
            "id": "smoke.impulse",
            "operation": "signal.impulse",
            "scipy_api": "numpy (constructed)",
            "versions": versions(),
            "params": {"length": n},
            "inputs": {},
            "expected": {"samples": arr(impulse)},
            "tolerance": {"rtol": 0.0, "atol": 0.0},
            "notes": "Unit impulse at index 0; base stimulus for convolution/filter fixtures.",
        }
    )
    write(
        {
            "id": "smoke.step",
            "operation": "signal.step",
            "scipy_api": "numpy (constructed)",
            "versions": versions(),
            "params": {"length": n},
            "inputs": {},
            "expected": {"samples": arr(step)},
            "tolerance": {"rtol": 0.0, "atol": 0.0},
        }
    )
    write(
        {
            "id": "smoke.sine_on_bin",
            "operation": "signal.sine",
            "scipy_api": "numpy.sin",
            "versions": versions(),
            "params": {
                "length": n,
                "sample_rate_hz": fs,
                "frequency_hz": 4.0 * fs / n,
            },
            "inputs": {},
            "expected": {"samples": arr(sine_on_bin)},
            "tolerance": {"rtol": 1e-15, "atol": 1e-15},
            "notes": "Tone exactly on an FFT bin for length 32.",
        }
    )
    write(
        {
            "id": "smoke.sine_between_bins",
            "operation": "signal.sine",
            "scipy_api": "numpy.sin",
            "versions": versions(),
            "params": {
                "length": n,
                "sample_rate_hz": fs,
                "frequency_hz": 12.3,
            },
            "inputs": {},
            "expected": {"samples": arr(sine_between)},
            "tolerance": {"rtol": 1e-15, "atol": 1e-15},
        }
    )

    x = impulse
    h_odd = np.array([0.25, 0.5, 0.25])
    h_even = np.array([0.1, 0.4, 0.4, 0.1])
    write(
        {
            "id": "smoke.convolve_full_causal_odd",
            "operation": "convolve.full",
            "scipy_api": "scipy.signal.convolve",
            "versions": versions(),
            "params": {"mode": "full"},
            "inputs": {"signal": arr(x), "kernel": arr(h_odd)},
            "expected": {"samples": arr(signal.convolve(x, h_odd, mode="full"))},
            "tolerance": {"rtol": 1e-12, "atol": 1e-14},
            "mapping": {
                "signal4s_kernel": {"zero_lag_index": 0, "constructor": "causal"},
                "signal4s_region": "Full",
            },
        }
    )
    write(
        {
            "id": "smoke.convolve_full_even_kernel_origin0",
            "operation": "convolve.full",
            "scipy_api": "scipy.signal.convolve",
            "versions": versions(),
            "params": {"mode": "full"},
            "inputs": {"signal": arr(x), "kernel": arr(h_even)},
            "expected": {"samples": arr(signal.convolve(x, h_even, mode="full"))},
            "tolerance": {"rtol": 1e-12, "atol": 1e-14},
            "mapping": {
                "signal4s_kernel": {"zero_lag_index": 0, "constructor": "at"},
                "signal4s_region": "Full",
                "notes": "Even-length kernel; origin must be explicit in signal4s.",
            },
        }
    )
    write(
        {
            "id": "smoke.convolve_valid_causal_odd",
            "operation": "convolve.valid",
            "scipy_api": "scipy.signal.convolve",
            "versions": versions(),
            "params": {"mode": "valid"},
            "inputs": {"signal": arr(x), "kernel": arr(h_odd)},
            "expected": {"samples": arr(signal.convolve(x, h_odd, mode="valid"))},
            "tolerance": {"rtol": 1e-12, "atol": 1e-14},
            "mapping": {
                "signal4s_kernel": {"zero_lag_index": 0, "constructor": "causal"},
                "signal4s_region": "Valid",
            },
        }
    )
    # SciPy mode='same' with odd kernel matches Input(Zero) for causal? For odd
    # length, SciPy centers on the middle tap. Map that to centeredOdd + Input(Zero).
    write(
        {
            "id": "smoke.convolve_same_odd_as_centered_input_zero",
            "operation": "convolve.input",
            "scipy_api": "scipy.signal.convolve",
            "versions": versions(),
            "params": {"mode": "same"},
            "inputs": {"signal": arr(x), "kernel": arr(h_odd)},
            "expected": {"samples": arr(signal.convolve(x, h_odd, mode="same"))},
            "tolerance": {"rtol": 1e-12, "atol": 1e-14},
            "mapping": {
                "signal4s_kernel": {
                    "zero_lag_index": 1,
                    "constructor": "centeredOdd",
                },
                "signal4s_region": "Input",
                "signal4s_boundary": "Zero",
                "notes": (
                    "SciPy same centers an odd kernel on the middle tap. "
                    "signal4s requires that origin explicitly via centeredOdd."
                ),
            },
        }
    )

    def correlate_signal4s(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, int]:
        """r[ℓ] = Σ_n a[n] b[n+ℓ]; returns values and first lag."""
        n = a.size
        m = b.size
        first = 1 - m
        last = n - 1
        out = np.empty(last - first + 1, dtype=np.float64)
        for i, lag in enumerate(range(first, last + 1)):
            start = max(0, -lag)
            end = min(n, m - lag)
            out[i] = (
                np.dot(a[start:end], b[start + lag : end + lag]) if end > start else 0.0
            )
        return out, first

    corr_vals, corr_first = correlate_signal4s(x, h_odd)
    write(
        {
            "id": "smoke.correlate_raw_impulse_odd",
            "operation": "correlate.raw",
            "scipy_api": "signal4s convention (see mapping); scipy.signal.correlate differs by lag sign",
            "versions": versions(),
            "params": {"normalization": "raw"},
            "inputs": {"x": arr(x), "y": arr(h_odd)},
            "expected": {
                "values": arr(corr_vals),
                "first_lag": corr_first,
            },
            "tolerance": {"rtol": 1e-12, "atol": 1e-14},
            "mapping": {
                "convention": "r_xy[lag] = sum_n x[n] * y[n+lag]",
                "scipy_note": (
                    "scipy.signal.correlate(x,y) uses R(k)=sum_n x[n+k] y[n], "
                    "i.e. signal4s lag = -scipy_lag under aligned indexing"
                ),
            },
        }
    )

    def write_fft_case(fid: str, x: np.ndarray, norm: str, kind: str) -> None:
        """kind: fft | rfft"""
        if kind == "fft":
            y = np.fft.fft(x, norm=norm)
            write(
                {
                    "id": fid,
                    "operation": "fft.forward",
                    "scipy_api": "numpy.fft.fft",
                    "versions": versions(),
                    "params": {"normalization": norm, "length": int(x.size)},
                    "inputs": {"real": arr(x.real), "imag": arr(np.zeros_like(x))},
                    "expected": {"real": arr(y.real), "imag": arr(y.imag)},
                    "tolerance": {"rtol": 1e-10, "atol": 1e-12},
                    "mapping": {
                        "signal4s_normalization": {
                            "backward": "Backward",
                            "forward": "Forward",
                            "ortho": "Orthonormal",
                        }[norm]
                    },
                }
            )
        else:
            y = np.fft.rfft(x, norm=norm)
            write(
                {
                    "id": fid,
                    "operation": "rfft.forward",
                    "scipy_api": "numpy.fft.rfft",
                    "versions": versions(),
                    "params": {"normalization": norm, "length": int(x.size)},
                    "inputs": {"signal": arr(x)},
                    "expected": {"real": arr(y.real), "imag": arr(y.imag)},
                    "tolerance": {"rtol": 1e-10, "atol": 1e-12},
                    "mapping": {
                        "signal4s_normalization": {
                            "backward": "Backward",
                            "forward": "Forward",
                            "ortho": "Orthonormal",
                        }[norm],
                        "layout": "NumPy rfft: n//2+1 complex bins",
                    },
                }
            )

    def write_ifft_case(fid: str, x: np.ndarray, norm: str) -> None:
        spectrum = np.fft.fft(x, norm=norm)
        y = np.fft.ifft(spectrum, norm=norm)
        write(
            {
                "id": fid,
                "operation": "fft.inverse",
                "scipy_api": "numpy.fft.ifft",
                "versions": versions(),
                "params": {"normalization": norm, "length": int(x.size)},
                "inputs": {"real": arr(spectrum.real), "imag": arr(spectrum.imag)},
                "expected": {"real": arr(y.real), "imag": arr(y.imag)},
                "tolerance": {"rtol": 1e-10, "atol": 1e-12},
                "mapping": {
                    "signal4s_normalization": {
                        "backward": "Backward",
                        "forward": "Forward",
                        "ortho": "Orthonormal",
                    }[norm]
                },
            }
        )

    def write_irfft_case(fid: str, x: np.ndarray, norm: str) -> None:
        spectrum = np.fft.rfft(x, norm=norm)
        y = np.fft.irfft(spectrum, n=x.size, norm=norm)
        write(
            {
                "id": fid,
                "operation": "rfft.inverse",
                "scipy_api": "numpy.fft.irfft",
                "versions": versions(),
                "params": {"normalization": norm, "length": int(x.size)},
                "inputs": {"real": arr(spectrum.real), "imag": arr(spectrum.imag)},
                "expected": {"signal": arr(y)},
                "tolerance": {"rtol": 1e-10, "atol": 1e-12},
                "mapping": {
                    "signal4s_normalization": {
                        "backward": "Backward",
                        "forward": "Forward",
                        "ortho": "Orthonormal",
                    }[norm],
                    "layout": "NumPy rfft: n//2+1 complex bins",
                },
            }
        )

    write_fft_case("smoke.fft_forward_sine_on_bin", sine_on_bin, "backward", "fft")
    write_fft_case("smoke.fft_forward_norm_ortho_16", sine_on_bin[:16], "ortho", "fft")
    write_fft_case(
        "smoke.fft_forward_norm_forward_16", sine_between[:16], "forward", "fft"
    )
    write_fft_case(
        "smoke.fft_forward_mixed_radix_12", sine_between[:12], "backward", "fft"
    )
    write_fft_case("smoke.fft_forward_prime_17", sine_between[:17], "backward", "fft")
    write_ifft_case("smoke.fft_inverse_backward_32", sine_between, "backward")
    write_ifft_case("smoke.fft_inverse_forward_16", sine_between[:16], "forward")
    write_fft_case("smoke.rfft_forward_sine_on_bin", sine_on_bin, "backward", "rfft")
    write_fft_case("smoke.rfft_forward_prime_17", sine_between[:17], "forward", "rfft")
    write_fft_case("smoke.rfft_forward_norm_ortho_32", sine_between, "ortho", "rfft")
    write_irfft_case("smoke.rfft_inverse_forward_17", sine_between[:17], "forward")
    write_irfft_case("smoke.rfft_inverse_ortho_32", sine_between, "ortho")

    # --- E4 filter fixtures ---
    fir_b = np.array([0.25, 0.5, 0.25], dtype=np.float64)
    fir_a = np.array([1.0], dtype=np.float64)
    x_f = sine_between
    y_fir = signal.lfilter(fir_b, fir_a, x_f)
    zi_fir = signal.lfilter_zi(fir_b, fir_a)
    y_fir_zi, zf_fir = signal.lfilter(fir_b, fir_a, x_f, zi=zi_fir * x_f[0])
    write(
        {
            "id": "smoke.lfilter_fir_causal",
            "operation": "lfilter.fir",
            "scipy_api": "scipy.signal.lfilter",
            "versions": versions(),
            "params": {},
            "inputs": {"signal": arr(x_f), "b": arr(fir_b), "a": arr(fir_a)},
            "expected": {"samples": arr(y_fir)},
            "tolerance": {"rtol": 1e-12, "atol": 1e-14},
            "mapping": {"signal4s": "Fir.causal(b).process(signal)"},
        }
    )
    write(
        {
            "id": "smoke.lfilter_fir_zi_zf",
            "operation": "lfilter.fir.zi",
            "scipy_api": "scipy.signal.lfilter / lfilter_zi",
            "versions": versions(),
            "params": {"zi_scale": "x[0]"},
            "inputs": {
                "signal": arr(x_f),
                "b": arr(fir_b),
                "a": arr(fir_a),
                "zi": arr(zi_fir * x_f[0]),
            },
            "expected": {"samples": arr(y_fir_zi), "zf": arr(zf_fir)},
            "tolerance": {"rtol": 1e-12, "atol": 1e-14},
            "mapping": {"signal4s": "FirRunner.restore(zi); process; snapshot ~= zf"},
        }
    )

    sos = signal.butter(4, 0.2, output="sos")
    y_sos = signal.sosfilt(sos, x_f)
    write(
        {
            "id": "smoke.sosfilt_butter4",
            "operation": "sosfilt",
            "scipy_api": "scipy.signal.sosfilt / butter",
            "versions": versions(),
            "params": {"order": 4, "wn": 0.2, "n_sections": int(sos.shape[0])},
            "inputs": {
                "signal": arr(x_f),
                "sos_flat": arr(sos.reshape(-1)),
            },
            "expected": {"samples": arr(y_sos)},
            "tolerance": {"rtol": 1e-10, "atol": 1e-12},
            "mapping": {"signal4s": "SecondOrderCascade.fromSosMatrix"},
        }
    )

    # Near-unit-circle pole stress: high-Q lowpass-ish biquad
    # poles near e^{±jθ} with radius 0.99
    r = 0.99
    theta = 0.15
    a1 = -2.0 * r * np.cos(theta)
    a2 = r * r
    b0, b1, b2 = 0.05, 0.1, 0.05
    sos_near = np.array([[b0, b1, b2, 1.0, a1, a2]], dtype=np.float64)
    y_near = signal.sosfilt(sos_near, x_f)
    write(
        {
            "id": "smoke.sosfilt_near_unit_circle",
            "operation": "sosfilt",
            "scipy_api": "scipy.signal.sosfilt",
            "versions": versions(),
            "params": {
                "pole_radius": r,
                "theta": theta,
                "n_sections": 1,
            },
            "inputs": {
                "signal": arr(x_f),
                "sos_flat": arr(sos_near.reshape(-1)),
            },
            "expected": {"samples": arr(y_near)},
            "tolerance": {"rtol": 1e-9, "atol": 1e-11},
            "notes": "Poles near the unit circle; numerical stress fixture.",
        }
    )

    y_ff = signal.filtfilt(fir_b, fir_a, x_f, padtype="odd")
    write(
        {
            "id": "smoke.filtfilt_fir_odd",
            "operation": "filtfilt",
            "scipy_api": "scipy.signal.filtfilt",
            "versions": versions(),
            "params": {"padtype": "odd", "method": "pad"},
            "inputs": {"signal": arr(x_f), "b": arr(fir_b), "a": arr(fir_a)},
            "expected": {"samples": arr(y_ff)},
            "tolerance": {"rtol": 1e-9, "atol": 1e-11},
            "mapping": {
                "signal4s": "ZeroPhase.filter(..., EdgeTreatment.OddPad)",
                "notes": "Matches SciPy padtype='odd' with default padlen.",
            },
        }
    )

    # --- E7 windows / STFT / Welch ---
    for name, n in (("hann", 16), ("hamming", 16), ("blackman", 16)):
        for fftbins, conv in ((True, "Periodic"), (False, "Symmetric")):
            w = signal.get_window(name, n, fftbins=fftbins)
            write(
                {
                    "id": f"smoke.window_{name}_{conv.lower()}_{n}",
                    "operation": "window",
                    "scipy_api": "scipy.signal.get_window",
                    "versions": versions(),
                    "params": {
                        "name": name,
                        "length": n,
                        "fftbins": fftbins,
                        "convention": conv,
                    },
                    "inputs": {},
                    "expected": {"taps": arr(w)},
                    "tolerance": {"rtol": 1e-12, "atol": 1e-14},
                    "mapping": {
                        "signal4s": f"WindowSpec.{name.capitalize()}(n, WindowConvention.{conv})",
                        "fftbins_true": "Periodic",
                        "fftbins_false": "Symmetric",
                    },
                }
            )

    k_sym = signal.windows.kaiser(16, beta=5.0, sym=True)
    k_per = signal.windows.kaiser(16, beta=5.0, sym=False)
    write(
        {
            "id": "smoke.window_kaiser_symmetric_16",
            "operation": "window",
            "scipy_api": "scipy.signal.windows.kaiser",
            "versions": versions(),
            "params": {
                "beta": 5.0,
                "length": 16,
                "sym": True,
                "convention": "Symmetric",
            },
            "inputs": {},
            "expected": {"taps": arr(k_sym)},
            "tolerance": {"rtol": 1e-7, "atol": 1e-9},
            "mapping": {
                "signal4s": "WindowSpec.Kaiser(16, 5.0, Symmetric)",
                "notes": "I0 via Abramowitz polynomial; looser tol than cosine windows.",
            },
        }
    )
    write(
        {
            "id": "smoke.window_kaiser_periodic_16",
            "operation": "window",
            "scipy_api": "scipy.signal.windows.kaiser",
            "versions": versions(),
            "params": {
                "beta": 5.0,
                "length": 16,
                "sym": False,
                "convention": "Periodic",
            },
            "inputs": {},
            "expected": {"taps": arr(k_per)},
            "tolerance": {"rtol": 1e-7, "atol": 1e-9},
            "mapping": {
                "signal4s": "WindowSpec.Kaiser(16, 5.0, Periodic)",
                "notes": "I0 via Abramowitz polynomial; looser tol than cosine windows.",
            },
        }
    )

    from scipy.signal import ShortTimeFFT

    x_stft = sine_between  # length 32
    # longer tone for STFT round-trip
    n_stft = 128
    fs_stft = 1000.0
    t_stft = np.arange(n_stft) / fs_stft
    x_stft = np.sin(2.0 * np.pi * 40.0 * t_stft)
    win_stft = signal.windows.hann(32, sym=False)
    hop = 16
    # phase_shift=None → raw windowed rFFT (SciPy default phase_shift=0 spins to window center).
    sft = ShortTimeFFT(
        win_stft, hop=hop, fs=fs_stft, fft_mode="onesided", phase_shift=None
    )
    Sx = sft.stft(x_stft)
    x_rt = sft.istft(Sx)[:n_stft]
    re = np.real(Sx).T.reshape(-1)  # frame-major: each frame's bins
    im = np.imag(Sx).T.reshape(-1)
    write(
        {
            "id": "smoke.stft_hann_periodic_hop16",
            "operation": "stft.analyze",
            "scipy_api": "scipy.signal.ShortTimeFFT",
            "versions": versions(),
            "params": {
                "sample_rate_hz": fs_stft,
                "frame_length": 32,
                "hop": hop,
                "nfft": 32,
                "window": "hann",
                "convention": "Periodic",
                "alignment": "Centered",
                "phase_shift": None,
                "n_frames": int(Sx.shape[1]),
                "n_bins": int(Sx.shape[0]),
            },
            "inputs": {"signal": arr(x_stft)},
            "expected": {
                "stft_real": arr(re),
                "stft_imag": arr(im),
                "roundtrip": arr(x_rt),
                "dual": arr(sft.dual_win),
            },
            "tolerance": {"rtol": 1e-9, "atol": 1e-11},
            "mapping": {
                "signal4s": "StftPlan.hannPeriodic(32, fs, hop=16)",
                "layout": "columns flattened in frame-major order (frame0 bins, frame1 bins, ...)",
                "notes": (
                    "Matches ShortTimeFFT(onesided, phase_shift=None). "
                    "SciPy's default phase_shift=0 applies a centering phase we do not."
                ),
            },
        }
    )

    # Welch / periodogram
    x_w = x_stft
    f_w, p_dens = signal.welch(
        x_w,
        fs=fs_stft,
        window=win_stft,
        nperseg=32,
        noverlap=16,
        nfft=32,
        detrend=False,
        return_onesided=True,
        scaling="density",
        average="mean",
    )
    write(
        {
            "id": "smoke.welch_density_nodetrend",
            "operation": "welch",
            "scipy_api": "scipy.signal.welch",
            "versions": versions(),
            "params": {
                "sample_rate_hz": fs_stft,
                "nperseg": 32,
                "noverlap": 16,
                "nfft": 32,
                "detrend": "none",
                "scaling": "density",
                "average": "mean",
                "sides": "onesided",
                "window": "hann_periodic",
            },
            "inputs": {"signal": arr(x_w)},
            "expected": {"frequencies": arr(f_w), "power": arr(p_dens)},
            "tolerance": {"rtol": 1e-9, "atol": 1e-12},
            "mapping": {
                "signal4s": "WelchPlan(..., detrend=None, scaling=Density)",
            },
        }
    )

    f_ws, p_spec = signal.welch(
        x_w,
        fs=fs_stft,
        window=win_stft,
        nperseg=32,
        noverlap=16,
        nfft=32,
        detrend="constant",
        return_onesided=True,
        scaling="spectrum",
        average="mean",
    )
    write(
        {
            "id": "smoke.welch_spectrum_detrend_mean",
            "operation": "welch",
            "scipy_api": "scipy.signal.welch",
            "versions": versions(),
            "params": {
                "sample_rate_hz": fs_stft,
                "nperseg": 32,
                "noverlap": 16,
                "nfft": 32,
                "detrend": "mean",
                "scaling": "spectrum",
                "average": "mean",
                "sides": "onesided",
                "window": "hann_periodic",
            },
            "inputs": {"signal": arr(x_w)},
            "expected": {"frequencies": arr(f_ws), "power": arr(p_spec)},
            "tolerance": {"rtol": 1e-9, "atol": 1e-12},
            "mapping": {
                "signal4s": "WelchPlan(..., detrend=Mean, scaling=Spectrum)",
            },
        }
    )

    f_p, p_p = signal.periodogram(
        x_w,
        fs=fs_stft,
        window="hann",
        nfft=128,
        detrend=False,
        return_onesided=True,
        scaling="spectrum",
    )
    write(
        {
            "id": "smoke.periodogram_spectrum_hann",
            "operation": "periodogram",
            "scipy_api": "scipy.signal.periodogram",
            "versions": versions(),
            "params": {
                "sample_rate_hz": fs_stft,
                "nfft": 128,
                "detrend": "none",
                "scaling": "spectrum",
                "sides": "onesided",
                "window": "hann",
            },
            "inputs": {"signal": arr(x_w)},
            "expected": {"frequencies": arr(f_p), "power": arr(p_p)},
            "tolerance": {"rtol": 1e-9, "atol": 1e-12},
            "mapping": {
                "signal4s": "Periodogram.hann(..., scaling=Spectrum, detrend=None)",
                "notes": "SciPy get_window('hann', len(x)) uses symmetric=False for spectral windows when fftbins default in get_window is True.",
            },
        }
    )

    # --- E8 filter design ---
    fs_d = 1000.0
    fir_lp = signal.firwin(11, 100.0, fs=fs_d, window="hamming", pass_zero=True)
    write(
        {
            "id": "smoke.firwin_lowpass_hamming_11",
            "operation": "firwin",
            "scipy_api": "scipy.signal.firwin",
            "versions": versions(),
            "params": {
                "numtaps": 11,
                "cutoff_hz": 100.0,
                "sample_rate_hz": fs_d,
                "window": "hamming",
                "pass_zero": True,
            },
            "inputs": {},
            "expected": {"taps": arr(fir_lp)},
            "tolerance": {"rtol": 1e-10, "atol": 1e-12},
            "mapping": {
                "signal4s": "FirDesign.lowPass(11, Frequency(100), SampleRate(1000), Hamming/Symmetric)",
            },
        }
    )

    z_b, p_b, k_b = signal.butter(4, 100.0, btype="low", output="zpk", fs=fs_d)
    sos_b = signal.butter(4, 100.0, btype="low", output="sos", fs=fs_d)
    ba_b = signal.butter(4, 100.0, btype="low", output="ba", fs=fs_d)
    w_f, h_f = signal.freqz(ba_b[0], ba_b[1], worN=16, fs=fs_d)
    _, h_s = signal.sosfreqz(sos_b, worN=16, fs=fs_d)
    t_i, y_i = signal.dimpulse((ba_b[0], ba_b[1], 1.0), n=16)
    t_s, y_s = signal.dstep((ba_b[0], ba_b[1], 1.0), n=16)
    write(
        {
            "id": "smoke.butter_lowpass4_100hz",
            "operation": "butter",
            "scipy_api": "scipy.signal.butter / freqz / sosfreqz / dimpulse / dstep",
            "versions": versions(),
            "params": {
                "order": 4,
                "cutoff_hz": 100.0,
                "sample_rate_hz": fs_d,
                "btype": "low",
            },
            "inputs": {},
            "expected": {
                "zpk_zeros_real": arr(np.real(z_b)),
                "zpk_zeros_imag": arr(np.imag(z_b)),
                "zpk_poles_real": arr(np.real(p_b)),
                "zpk_poles_imag": arr(np.imag(p_b)),
                "zpk_gain": float(k_b),
                "sos_flat": arr(sos_b.reshape(-1)),
                "freqz_mag": arr(np.abs(h_f)),
                "sosfreqz_mag": arr(np.abs(h_s)),
                "impulse": arr(y_i[0].ravel()),
                "step": arr(y_s[0].ravel()),
            },
            "tolerance": {"rtol": 1e-9, "atol": 1e-11},
            "mapping": {
                "signal4s": "Butterworth.lowPass(4, Frequency(100), SampleRate(1000))",
                "notes": (
                    "ZPK must match SciPy. SOS section order/pairing may differ; "
                    "compare frequency and impulse responses rather than row-wise SOS."
                ),
            },
        }
    )

    n_k, beta_k = signal.kaiserord(40, 50.0 / (fs_d / 2.0))
    write(
        {
            "id": "smoke.kaiserord_40db_50hz",
            "operation": "kaiserord",
            "scipy_api": "scipy.signal.kaiserord",
            "versions": versions(),
            "params": {
                "ripple_db": 40.0,
                "width_hz": 50.0,
                "sample_rate_hz": fs_d,
            },
            "inputs": {},
            "expected": {"numtaps": int(n_k), "beta": float(beta_k)},
            "tolerance": {"rtol": 0.0, "atol": 1e-12},
            "mapping": {"signal4s": "FirDesign.kaiserOrder(40, 50, SampleRate(1000))"},
        }
    )

    # --- E9 multirate ---
    h_u = np.array([0.25, 0.5, 0.25], dtype=np.float64)
    x_u = np.arange(1, 9, dtype=np.float64)
    for up, down in ((1, 1), (2, 1), (1, 2), (2, 3), (3, 2)):
        y_u = signal.upfirdn(h_u, x_u, up=up, down=down)
        write(
            {
                "id": f"smoke.upfirdn_up{up}_down{down}",
                "operation": "upfirdn",
                "scipy_api": "scipy.signal.upfirdn",
                "versions": versions(),
                "params": {"up": up, "down": down},
                "inputs": {"h": arr(h_u), "signal": arr(x_u)},
                "expected": {"samples": arr(y_u)},
                "tolerance": {"rtol": 1e-12, "atol": 1e-14},
                "mapping": {
                    "signal4s": f"Upfirdn(h, x, up={up}, down={down})",
                    "notes": (
                        "Zero-insert between samples only (not after last), "
                        "full FIR, keep every down-th sample from index 0. "
                        "Streaming PolyphaseResampler.consume+flush matches this batch."
                    ),
                },
            }
        )

    x_r = np.sin(2.0 * np.pi * np.arange(32) * 3.0 / 32.0)
    for up, down in ((3, 2), (2, 3), (4, 2)):
        y_r = signal.resample_poly(x_r, up, down, window=("kaiser", 5.0))
        write(
            {
                "id": f"smoke.resample_poly_up{up}_down{down}",
                "operation": "resample_poly",
                "scipy_api": "scipy.signal.resample_poly",
                "versions": versions(),
                "params": {
                    "up": up,
                    "down": down,
                    "window": ["kaiser", 5.0],
                    "padtype": "constant",
                },
                "inputs": {"signal": arr(x_r)},
                "expected": {"samples": arr(y_r)},
                "tolerance": {"rtol": 1e-8, "atol": 1e-10},
                "mapping": {
                    "signal4s": "ResamplePoly(signal, up, down)",
                    "notes": "Includes SciPy filter *up gain, pre-pad centering, and output trim.",
                },
            }
        )

    # --- E10 end-to-end: generate → butter SOS → Welch ---
    fs_e = 1000.0
    n_e = 256
    t_e = np.arange(n_e) / fs_e
    x_e = np.sin(2.0 * np.pi * 40.0 * t_e) + 0.25 * np.sin(2.0 * np.pi * 120.0 * t_e)
    sos_e = signal.butter(4, 80.0, btype="low", output="sos", fs=fs_e)
    y_e = signal.sosfilt(sos_e, x_e)
    f_e, p_e = signal.welch(
        y_e,
        fs=fs_e,
        window="hann",
        nperseg=64,
        noverlap=32,
        nfft=64,
        detrend="constant",
        return_onesided=True,
        scaling="density",
        average="mean",
    )
    write(
        {
            "id": "smoke.e2e_butter_welch",
            "operation": "e2e.filter_welch",
            "scipy_api": "butter(sos) → sosfilt → welch",
            "versions": versions(),
            "params": {
                "sample_rate_hz": fs_e,
                "butter_order": 4,
                "cutoff_hz": 80.0,
                "welch_nperseg": 64,
                "welch_noverlap": 32,
                "detrend": "mean",
                "scaling": "density",
            },
            "inputs": {
                "signal": arr(x_e),
                "sos_flat": arr(sos_e.reshape(-1)),
            },
            "expected": {
                "filtered": arr(y_e),
                "frequencies": arr(f_e),
                "power": arr(p_e),
            },
            "tolerance": {"rtol": 1e-8, "atol": 1e-10},
            "mapping": {
                "signal4s": (
                    "Butterworth.lowPass → sos.process → WelchPlan "
                    "(or sos_flat via SecondOrderCascade.fromSosMatrix for exact SOS match)"
                ),
                "notes": "End-to-end workflow fixture for the 1.0 audit.",
            },
        }
    )


if __name__ == "__main__":
    main()
