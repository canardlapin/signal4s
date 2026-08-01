#!/usr/bin/env python3
"""Time SciPy/NumPy peers for vs-signal4s receipts. Not invoked by CI.

Shapes are a *grid* over sizes and configurations—not a shortlist tuned to a
single optimization. Median of several trials after warm-up.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import scipy
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "benchmarks" / "receipts" / "vs-scipy-python.json"

WARMUP = 50
ITERS = 100
TRIALS = 7


def next_pow2(n: int) -> int:
    return 1 << (n - 1).bit_length()


def time_ms(
    fn, warmup: int = WARMUP, iters: int = ITERS, trials: int = TRIALS
) -> float:
    for _ in range(warmup):
        fn()
    samples: list[float] = []
    for _ in range(trials):
        t0 = time.perf_counter()
        sink = 0.0
        for _ in range(iters):
            y = fn()
            y0 = np.asarray(y).ravel()[0]
            sink += float(np.real(y0))
        _ = sink
        samples.append((time.perf_counter() - t0) * 1e3 / iters)
    samples.sort()
    return samples[len(samples) // 2]


def shapes() -> list[dict]:
    """Build one full receipt. Caller may discard a settle pass."""
    rows: list[dict] = []

    for n, m in (
        (128, 9),
        (256, 16),
        (512, 33),
        (1024, 48),
        (2048, 64),
        (4096, 80),
        (8192, 128),
        (3000, 50),
    ):
        x = np.sin(0.02 * np.arange(n, dtype=np.float64))
        h = 1.0 / (np.arange(m, dtype=np.float64) + 1.0)

        def run_direct(x=x, h=h):
            return signal.convolve(x, h, mode="full", method="direct")

        rows.append(
            {
                "id": f"convolve_direct_full_n{n}_m{m}",
                "op": "convolve_direct_full",
                "n": n,
                "m": m,
                "ms": time_ms(run_direct),
            }
        )

        def run_fft(x=x, h=h):
            return signal.convolve(x, h, mode="full", method="fft")

        rows.append(
            {
                "id": f"convolve_fft_full_n{n}_m{m}",
                "op": "convolve_fft_full",
                "n": n,
                "m": m,
                "ms": time_ms(run_fft),
            }
        )

        nfft = next_pow2(n + m - 1)
        kernel_spectrum = np.fft.rfft(h, n=nfft)

        def run_fft_planned(
            x=x, kernel_spectrum=kernel_spectrum, nfft=nfft, out_len=n + m - 1
        ):
            spectrum = np.fft.rfft(x, n=nfft)
            return np.fft.irfft(spectrum * kernel_spectrum, n=nfft)[:out_len]

        rows.append(
            {
                "id": f"convolve_fft_planned_full_n{n}_m{m}",
                "op": "convolve_fft_planned_full",
                "n": n,
                "m": m,
                "nfft": nfft,
                "ms": time_ms(run_fft_planned),
            }
        )

    for n, m in ((1024, 48), (4096, 80), (8192, 128)):
        x = np.sin(0.02 * np.arange(n, dtype=np.float64))
        h = 1.0 / (np.arange(m, dtype=np.float64) + 1.0)

        def run_fir_causal(x=x, h=h):
            return signal.lfilter(h, [1.0], x)

        rows.append(
            {
                "id": f"fir_causal_n{n}_m{m}",
                "op": "fir_causal",
                "n": n,
                "m": m,
                "ms": time_ms(run_fir_causal),
            }
        )

        def run_valid_direct(x=x, h=h):
            return signal.convolve(x, h, mode="valid", method="direct")

        rows.append(
            {
                "id": f"convolve_valid_direct_n{n}_m{m}",
                "op": "convolve_valid_direct",
                "n": n,
                "m": m,
                "ms": time_ms(run_valid_direct),
            }
        )

        kernel_circular = np.zeros(n, dtype=np.float64)
        kernel_circular[:m] = h

        def run_circular_fft(x=x, kernel_circular=kernel_circular, n=n):
            return np.fft.irfft(np.fft.rfft(x) * np.fft.rfft(kernel_circular), n=n)

        rows.append(
            {
                "id": f"circular_fft_n{n}_m{m}",
                "op": "circular_fft",
                "n": n,
                "m": m,
                "ms": time_ms(run_circular_fft),
            }
        )

    for n in (64, 256, 1024, 2048, 4096, 8192):
        x = np.sin(0.02 * np.arange(n, dtype=np.float64)) + 0j

        def run(x=x):
            return np.fft.fft(x)

        rows.append(
            {"id": f"fft_pow2_n{n}", "op": "fft_pow2", "n": n, "ms": time_ms(run)}
        )

        xr = np.sin(0.02 * np.arange(n, dtype=np.float64))

        def run_real(xr=xr):
            return np.fft.rfft(xr)

        rows.append(
            {
                "id": f"rfft_pow2_n{n}",
                "op": "rfft_pow2",
                "n": n,
                "ms": time_ms(run_real),
            }
        )

    for n in (12, 60, 2160, 8640):
        x = np.sin(0.02 * np.arange(n, dtype=np.float64)) + 1j * np.cos(
            0.017 * np.arange(n, dtype=np.float64)
        )

        def run_smooth(x=x):
            return np.fft.fft(x)

        rows.append(
            {
                "id": f"fft_smooth_n{n}",
                "op": "fft_smooth",
                "n": n,
                "ms": time_ms(run_smooth),
            }
        )

    sos = signal.butter(4, 0.2, output="sos")
    for n in (512, 2048, 8192, 32768, 65536):
        x = np.sin(0.02 * np.arange(n, dtype=np.float64))

        def run(x=x, sos=sos):
            return signal.sosfilt(sos, x)

        rows.append(
            {
                "id": f"sosfilt_n{n}_sec{sos.shape[0]}",
                "op": "sosfilt",
                "n": n,
                "sections": int(sos.shape[0]),
                "ms": time_ms(run),
            }
        )

    for n, m, up, down in (
        (512, 15, 2, 3),
        (2048, 31, 3, 2),
        (4096, 47, 4, 3),
        (8192, 63, 5, 3),
        (16384, 63, 2, 1),
        (8192, 63, 1, 2),
    ):
        x = np.sin(0.02 * np.arange(n, dtype=np.float64))
        h = np.hamming(m).astype(np.float64)

        def run(x=x, h=h, up=up, down=down):
            return signal.upfirdn(h, x, up=up, down=down)

        rows.append(
            {
                "id": f"upfirdn_n{n}_m{m}_up{up}_down{down}",
                "op": "upfirdn",
                "n": n,
                "m": m,
                "up": up,
                "down": down,
                "ms": time_ms(run),
            }
        )

    return rows


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    # First full pass settles caches; receipt uses the second pass only.
    _ = shapes()
    rows = shapes()
    payload = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "versions": {"scipy": scipy.__version__, "numpy": np.__version__},
        "warmup": WARMUP,
        "iters": ITERS,
        "trials": TRIALS,
        "rows": rows,
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"wrote {OUT}")
    for r in payload["rows"]:
        print(f"{r['id']}: {r['ms']:.4f} ms")


if __name__ == "__main__":
    main()
