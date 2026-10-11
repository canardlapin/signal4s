#!/usr/bin/env python3
"""Freeze one-frame STFT answers from pinned public SciPy, outside JVM tests."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import scipy
from scipy.signal import ShortTimeFFT, detrend as remove_trend
from scipy.signal.windows import hann

ROOT = Path(__file__).resolve().parents[2]
assert sys.flags.isolated, 'use isolated pinned reference interpreter'
assert scipy.__version__ == '1.18.1' and np.__version__ == '2.3.5'


def main():
    definitions = [
        ('impulse_phase', [0., 1., 0., 0.], np.ones(4), 4, 'None', 'Raw'),
        ('hann_analysis_only', [0., 0., 1., 0.], hann(4, sym=False), 4, 'None', 'Raw'),
        ('odd_prime_linear', [1., 3., -2., 5., 2.], np.ones(5), 7, 'Linear', 'Density'),
        ('prime_symmetric', [1., -2., 3., -4., 5., -6., 7.], hann(7, sym=True), 11, 'Mean', 'Raw'),
        ('padded_periodic', [1., 2., -3., 4., 5., -6., 7., 8.], hann(8, sym=False), 16, 'Mean', 'Magnitude'),
        ('signed_window', [1., 2., 3., 4.], [-1., -1., -1., -1.], 4, 'None', 'Magnitude'),
        ('zero_coherent', [1., 2., 3., 4.], [1., -1., 1., -1.], 5, 'None', 'Density'),
        ('tiny_window', [1., -2., 3., 4.], np.ones(4) * 1e-150, 8, 'None', 'Density'),
        ('large_window', [1., -2., 3., 4.], np.ones(4) * 1e150, 8, 'None', 'Magnitude'),
        ('large_input', np.array([1., -2., 3., 4.]) * 1e150, np.ones(4), 7, 'None', 'Density'),
        ('singleton', [2.], [1.], 1, 'None', 'Density'),
        ('tone_no_doubling', np.cos(2*np.pi*np.arange(8)/8), np.ones(8), 8, 'None', 'Magnitude'),
    ]
    cases = []
    for name, values, window, nfft, trend, scaling in definitions:
        x, w = np.asarray(values, float), np.asarray(window, float)
        mode = {'Raw': None, 'Magnitude': 'magnitude', 'Density': 'psd'}[scaling]
        transform = ShortTimeFFT(w, hop=len(w), fs=16., fft_mode='onesided', mfft=nfft,
                                 scale_to=mode, phase_shift=None)
        detr = {'None': None, 'Mean': 'constant', 'Linear': 'linear'}[trend]
        actual = transform.stft_detrend(x, detr, p0=0, p1=1, k_offset=len(w)//2)[:, 0]
        y = x if detr is None else remove_trend(x, type=detr)
        divisor = 1. if mode is None else abs(w.sum()) if mode == 'magnitude' else np.sqrt(16. * np.sum(w*w))
        independent = np.fft.rfft(y*w, n=nfft) / divisor
        np.testing.assert_allclose(actual, independent, rtol=1e-12, atol=1e-12*max(np.max(np.abs(independent)), 1e-300))
        cases.append({'id': name, 'input': x.tolist(), 'window': w.tolist(), 'nfft': nfft,
                      'sample_rate': 16., 'detrend': trend, 'scaling': scaling,
                      'real': actual.real.tolist(), 'imaginary': actual.imag.tolist(),
                      'frequency': transform.f.tolist()})
    report = {'format': 'signal4s-bounded-stft-goldens-v1', 'reference': {'scipy': scipy.__version__,
              'numpy': np.__version__, 'api': 'public ShortTimeFFT.stft_detrend; full frame, phase_shift=None, onesided'},
              'relative_tolerance': 1e-10, 'absolute_floor': 1e-300, 'cases': cases}
    dest = ROOT / 'fixtures/data/bounded-stft-v1.json'
    assert not dest.exists(), 'do not silently regenerate goldens'
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    def vector(values):
        return 'Vector(' + ','.join(repr(float(v)) for v in values) + ')'
    lines = ['package signal4s.fft', '', 'private[fft] object BoundedStftGoldens:',
             '  final case class Case(id:String,input:Vector[Double],window:Vector[Double],nfft:Int,',
             '      detrend:Detrend,scaling:StftScaling,real:Vector[Double],imaginary:Vector[Double])',
             '  val cases: Vector[Case] = Vector(']
    for index, c in enumerate(cases):
        lines.append('    Case(' + json.dumps(c['id']) + ',' + vector(c['input']) + ',' + vector(c['window']) +
                     ',' + str(c['nfft']) + ',Detrend.' + c['detrend'] + ',StftScaling.' + c['scaling'] +
                     ',' + vector(c['real']) + ',' + vector(c['imaginary']) + ')' + (',' if index+1 < len(cases) else ''))
    lines.append('  )')
    scala = ROOT / 'modules/fft/shared/src/test/scala/signal4s/fft/BoundedStftGoldens.scala'
    assert not scala.exists()
    scala.write_text('\n'.join(lines) + '\n')
    print('frozen', len(cases), 'public SciPy cases', hashlib.sha256(dest.read_bytes()).hexdigest())


if __name__ == '__main__':
    main()
