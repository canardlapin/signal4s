# E6 Auto cost model calibration

generated: 2026-08-01T20:39:04.812162Z

Planned wall-clock (ms/call) for Full convolution; warm runs.

| N | M | Direct | FFT | OLA(B) | Auto picks |
|---|---|--------|-----|--------|------------|
| 64 | 8 | 0.0145 | 0.0419 | 0.0240 (B=32) | Direct |
| 256 | 16 | 0.0169 | 0.0268 | 0.0207 (B=64) | Direct |
| 1024 | 32 | 0.0276 | 0.0484 | 0.0470 (B=128) | Direct |
| 4096 | 64 | 0.0686 | 0.1799 | 0.1876 (B=256) | Direct |
| 8192 | 9 | 0.1174 | 0.2440 | 0.1973 (B=64) | Direct |
| 1024 | 512 | 0.0399 | 0.0248 | 0.0269 (B=1024) | Fft |
| 8192 | 128 | 0.0986 | 0.2062 | 0.1201 (B=512) | Direct |
| 65536 | 64 | 0.4047 | 1.9065 | 0.7645 (B=256) | Direct |
| 16384 | 2048 | 8.0070 | 0.2926 | 0.2718 (B=8192) | OverlapAdd(8192) |
| 16384 | 4096 | 16.0602 | 0.2967 | 0.2919 (B=16384) | Fft |

## Model

AutoCostModel uses relative costs:
- Direct ≈ 1.0 · N · M
- FFT ≈ 12.0 · L · log2(L), L = nextPow2(N+M-1)
- OLA ≈ 12.0 · ⌈N/B⌉ · L_b · log2(L_b), L_b = nextPow2(B+M-1)

Coefficients chosen so ordering tracks this JVM portable matrix.
