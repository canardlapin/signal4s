# Performance vs SciPy

generated: 2026-08-01T16:20:04.081040Z

See [`docs/PERFORMANCE.md`](../../docs/PERFORMANCE.md) for fairness rules.

Trials: median of 7 × 100 iters after 50 warm-up.

_Python receipt missing at `benchmarks/receipts/vs-scipy-python.json`. Run `python3 fixtures/generate/bench_vs_scipy.py` first._
- `convolve_direct_full_n128_m9`: scala=0.0017 ms
- `convolve_fft_full_n128_m9`: scala=0.0016 ms
- `convolve_direct_full_n256_m16`: scala=0.0277 ms
- `convolve_fft_full_n256_m16`: scala=0.0033 ms
- `convolve_direct_full_n512_m33`: scala=0.1078 ms
- `convolve_fft_full_n512_m33`: scala=0.0067 ms
- `convolve_direct_full_n1024_m48`: scala=0.3218 ms
- `convolve_fft_full_n1024_m48`: scala=0.0148 ms
- `convolve_direct_full_n2048_m64`: scala=0.8164 ms
- `convolve_fft_full_n2048_m64`: scala=0.0296 ms
- `convolve_direct_full_n4096_m80`: scala=2.0871 ms
- `convolve_fft_full_n4096_m80`: scala=0.0739 ms
- `convolve_direct_full_n8192_m128`: scala=6.7346 ms
- `convolve_fft_full_n8192_m128`: scala=0.1536 ms
- `convolve_direct_full_n3000_m50`: scala=0.9687 ms
- `convolve_fft_full_n3000_m50`: scala=0.0321 ms
- `fft_pow2_n64`: scala=0.0004 ms
- `fft_pow2_n256`: scala=0.0013 ms
- `fft_pow2_n1024`: scala=0.0056 ms
- `fft_pow2_n2048`: scala=0.0115 ms
- `fft_pow2_n4096`: scala=0.0278 ms
- `fft_pow2_n8192`: scala=0.0608 ms
- `sosfilt_n512_sec2`: scala=0.0019 ms
- `sosfilt_n2048_sec2`: scala=0.0079 ms
- `sosfilt_n8192_sec2`: scala=0.0308 ms
- `sosfilt_n32768_sec2`: scala=0.1166 ms
- `sosfilt_n65536_sec2`: scala=0.2312 ms
- `upfirdn_n512_m15_up2_down3`: scala=0.0035 ms
- `upfirdn_n2048_m31_up3_down2`: scala=0.0262 ms
- `upfirdn_n4096_m47_up4_down3`: scala=0.0544 ms
- `upfirdn_n8192_m63_up5_down3`: scala=0.1408 ms
- `upfirdn_n16384_m63_up2_down1`: scala=0.3173 ms
- `upfirdn_n8192_m63_up1_down2`: scala=0.2085 ms
