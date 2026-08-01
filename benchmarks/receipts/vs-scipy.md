# Performance vs SciPy

generated: 2026-08-01T22:17:25.642886Z

See [`docs/PERFORMANCE.md`](../../docs/PERFORMANCE.md) for fairness rules.

Trials: median of 7 × 100 iters after 50 warm-up.

| id | scipy ms | signal4s ms | ratio (s4s/scipy) |
|---|---:|---:|---:|
| `convolve_direct_full_n128_m9` | 0.0025 | 0.0017 | 0.70× |
| `convolve_fft_full_n128_m9` | 0.0205 | 0.0029 | 0.14× |
| `convolve_fft_planned_full_n128_m9` | 0.0070 | 0.0015 | 0.21× |
| `convolve_direct_full_n256_m16` | 0.0036 | 0.0006 | 0.16× |
| `convolve_fft_full_n256_m16` | 0.0228 | 0.0051 | 0.22× |
| `convolve_fft_planned_full_n256_m16` | 0.0088 | 0.0030 | 0.34× |
| `convolve_direct_full_n512_m33` | 0.0052 | 0.0015 | 0.28× |
| `convolve_fft_full_n512_m33` | 0.0254 | 0.0109 | 0.43× |
| `convolve_fft_planned_full_n512_m33` | 0.0115 | 0.0060 | 0.52× |
| `convolve_direct_full_n1024_m48` | 0.0090 | 0.0034 | 0.38× |
| `convolve_fft_full_n1024_m48` | 0.0306 | 0.0220 | 0.72× |
| `convolve_fft_planned_full_n1024_m48` | 0.0184 | 0.0122 | 0.66× |
| `convolve_direct_full_n2048_m64` | 0.0199 | 0.0085 | 0.43× |
| `convolve_fft_full_n2048_m64` | 0.0413 | 0.0448 | 1.08× |
| `convolve_fft_planned_full_n2048_m64` | 0.0296 | 0.0274 | 0.93× |
| `convolve_direct_full_n4096_m80` | 0.0392 | 0.0196 | 0.50× |
| `convolve_fft_full_n4096_m80` | 0.0694 | 0.0961 | 1.38× |
| `convolve_fft_planned_full_n4096_m80` | 0.0765 | 0.0610 | 0.80× |
| `convolve_direct_full_n8192_m128` | 0.1131 | 0.0645 | 0.57× |
| `convolve_fft_full_n8192_m128` | 0.1128 | 0.2053 | 1.82× |
| `convolve_fft_planned_full_n8192_m128` | 0.1281 | 0.1274 | 0.99× |
| `convolve_direct_full_n3000_m50` | 0.0227 | 0.0101 | 0.44× |
| `convolve_fft_full_n3000_m50` | 0.0484 | 0.0452 | 0.93× |
| `convolve_fft_planned_full_n3000_m50` | 0.0295 | 0.0281 | 0.95× |
| `fir_causal_n1024_m48` | 0.0190 | 0.0033 | 0.18× |
| `convolve_valid_direct_n1024_m48` | 0.0092 | 0.0036 | 0.39× |
| `circular_fft_n1024_m48` | 0.0160 | 0.0103 | 0.64× |
| `fir_causal_n4096_m80` | 0.0517 | 0.0190 | 0.37× |
| `convolve_valid_direct_n4096_m80` | 0.0391 | 0.0192 | 0.49× |
| `circular_fft_n4096_m80` | 0.0438 | 0.0430 | 0.98× |
| `fir_causal_n8192_m128` | 0.1286 | 0.0537 | 0.42× |
| `convolve_valid_direct_n8192_m128` | 0.1132 | 0.0535 | 0.47× |
| `circular_fft_n8192_m128` | 0.0970 | 0.0942 | 0.97× |
| `fft_pow2_n64` | 0.0029 | 0.0004 | 0.14× |
| `rfft_pow2_n64` | 0.0028 | 0.0002 | 0.06× |
| `fft_pow2_n256` | 0.0037 | 0.0012 | 0.34× |
| `rfft_pow2_n256` | 0.0032 | 0.0007 | 0.21× |
| `fft_pow2_n1024` | 0.0070 | 0.0053 | 0.76× |
| `rfft_pow2_n1024` | 0.0052 | 0.0028 | 0.53× |
| `fft_pow2_n2048` | 0.0129 | 0.0117 | 0.91× |
| `rfft_pow2_n2048` | 0.0086 | 0.0056 | 0.66× |
| `fft_pow2_n4096` | 0.0245 | 0.0277 | 1.13× |
| `rfft_pow2_n4096` | 0.0134 | 0.0156 | 1.16× |
| `fft_pow2_n8192` | 0.0586 | 0.0626 | 1.07× |
| `rfft_pow2_n8192` | 0.0323 | 0.0384 | 1.19× |
| `fft_smooth_n12` | 0.0029 | 0.0001 | 0.03× |
| `fft_smooth_n60` | 0.0029 | 0.0003 | 0.12× |
| `fft_smooth_n2160` | 0.0145 | 0.0162 | 1.12× |
| `fft_smooth_n8640` | 0.0564 | 0.0821 | 1.46× |
| `sosfilt_n512_sec2` | 0.0110 | 0.0032 | 0.29× |
| `sosfilt_n2048_sec2` | 0.0179 | 0.0126 | 0.70× |
| `sosfilt_n8192_sec2` | 0.0451 | 0.0291 | 0.64× |
| `sosfilt_n32768_sec2` | 0.1546 | 0.1144 | 0.74× |
| `sosfilt_n65536_sec2` | 0.3203 | 0.2294 | 0.72× |
| `upfirdn_n512_m15_up2_down3` | 0.0048 | 0.0023 | 0.48× |
| `upfirdn_n2048_m31_up3_down2` | 0.0184 | 0.0162 | 0.88× |
| `upfirdn_n4096_m47_up4_down3` | 0.0321 | 0.0323 | 1.01× |
| `upfirdn_n8192_m63_up5_down3` | 0.0793 | 0.0801 | 1.01× |
| `upfirdn_n16384_m63_up2_down1` | 0.4662 | 0.3121 | 0.67× |
| `upfirdn_n8192_m63_up1_down2` | 0.1465 | 0.2075 | 1.42× |
