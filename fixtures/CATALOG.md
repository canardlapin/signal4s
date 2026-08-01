# Fixture coverage catalog

Maps 1.0 required numerical features to fixture ids. Status **present + tested**
means a JVM suite loads the committed JSON (no Python at test time).

| Feature | Fixture ids | SciPy API | Epic | Status |
|---|---|---|---|---|
| Impulse / step / sine stimuli | `smoke.impulse`, `smoke.step`, `smoke.sine_on_bin`, `smoke.sine_between_bins` | numpy constructed | E0 | present |
| Convolve full (causal odd) | `smoke.convolve_full_causal_odd` | `scipy.signal.convolve` | E2 | present + tested |
| Convolve full (even kernel, explicit origin) | `smoke.convolve_full_even_kernel_origin0` | `scipy.signal.convolve` | E2 | present + tested |
| Convolve valid | `smoke.convolve_valid_causal_odd` | `scipy.signal.convolve` | E2 | present + tested |
| Convolve same → Input(Zero) via centeredOdd | `smoke.convolve_same_odd_as_centered_input_zero` | `scipy.signal.convolve` | E2 | present + tested |
| Correlate + lag axis (signal4s convention) | `smoke.correlate_raw_impulse_odd` | explicit convention (SciPy lag sign differs) | E2 | present + tested |
| FIR lfilter | `smoke.lfilter_fir_causal` | `scipy.signal.lfilter` | E4 | present + tested |
| FIR zi/zf | `smoke.lfilter_fir_zi_zf` | `lfilter` / `lfilter_zi` | E4 | present + tested |
| SOS butter | `smoke.sosfilt_butter4` | `sosfilt` / `butter` | E4 | present + tested |
| SOS near unit circle | `smoke.sosfilt_near_unit_circle` | `sosfilt` | E4 | present + tested |
| filtfilt odd pad | `smoke.filtfilt_fir_odd` | `filtfilt(padtype='odd')` | E4 | present + tested |
| FFT forward (backward) | `smoke.fft_forward_sine_on_bin` | `numpy.fft.fft` | E5 | present + tested |
| FFT forward (forward / ortho) | `smoke.fft_forward_norm_forward_16`, `smoke.fft_forward_norm_ortho_16` | `numpy.fft.fft(norm=...)` | E5 | present + tested |
| FFT inverse (backward / forward) | `smoke.fft_inverse_backward_32`, `smoke.fft_inverse_forward_16` | `numpy.fft.ifft(norm=...)` | E5 | present + tested |
| FFT mixed-radix 12 | `smoke.fft_forward_mixed_radix_12` | `numpy.fft.fft` | E5 | present + tested |
| FFT Bluestein prime 17 | `smoke.fft_forward_prime_17` | `numpy.fft.fft` | E5 | present + tested |
| rfft backward | `smoke.rfft_forward_sine_on_bin` | `numpy.fft.rfft` | E5 | present + tested |
| rfft forward (forward / ortho) | `smoke.rfft_forward_prime_17`, `smoke.rfft_forward_norm_ortho_32` | `numpy.fft.rfft(norm=...)` | E5 | present + tested |
| rfft inverse (forward / ortho) | `smoke.rfft_inverse_forward_17`, `smoke.rfft_inverse_ortho_32` | `numpy.fft.irfft(norm=...)` | E5 | present + tested |
| Windows (Hann/Hamming/Blackman/Kaiser × Periodic/Symmetric) | `smoke.window_*` | `get_window` / `windows.kaiser` | E7 | present + tested |
| STFT + dual + round-trip | `smoke.stft_hann_periodic_hop16` | `ShortTimeFFT(phase_shift=None)` | E7 | present + tested |
| Welch density / spectrum | `smoke.welch_density_nodetrend`, `smoke.welch_spectrum_detrend_mean` | `scipy.signal.welch` | E7 | present + tested |
| Periodogram spectrum | `smoke.periodogram_spectrum_hann` | `scipy.signal.periodogram` | E7 | present + tested |
| firwin lowpass Hamming | `smoke.firwin_lowpass_hamming_11` | `scipy.signal.firwin` | E8 | present + tested |
| Kaiser order | `smoke.kaiserord_40db_50hz` | `scipy.signal.kaiserord` | E8 | present + tested |
| Butterworth LP ZPK/freqz/impulse/step | `smoke.butter_lowpass4_100hz` | `butter` / `freqz` / `sosfreqz` / `dimpulse` / `dstep` | E8 | present + tested |
| upfirdn (several up/down) | `smoke.upfirdn_up*_down*` | `scipy.signal.upfirdn` | E9 | present + tested |
| resample_poly | `smoke.resample_poly_up*_down*` | `scipy.signal.resample_poly` | E9 | present + tested |
| E2E butter SOS → Welch | `smoke.e2e_butter_welch` | `butter`/`sosfilt`/`welch` | E10 | present + tested |

**Coverage:** 100% of required §6 numerical features for 1.0.

Pinned generator versions: see `fixtures/generate/requirements.txt`
(`numpy==2.2.6`, `scipy==1.15.3`).

Intentional semantic differences: [`docs/SEMANTICS.md`](../docs/SEMANTICS.md).
