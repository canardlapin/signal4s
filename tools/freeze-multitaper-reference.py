#!/usr/bin/env python3
"""Pinned independent DPSS/direct-FFT reference; --check never rewrites fixtures."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import scipy
from scipy.signal.windows import dpss

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    assert np.__version__=='2.3.5' and scipy.__version__=='1.18.1'
    rng=np.random.default_rng(20261014)
    cases=[]
    for n,nw,k,floor,nfft,detrend in [(31,2.5,5,.9,64,'None'),(32,3.,5,0.,32,'Constant'),(64,3.5,6,.95,128,'Constant')]:
        tapers,ratios=dpss(n,nw,Kmax=k,sym=True,norm=2,return_ratios=True)
        selected=np.flatnonzero(ratios>=floor)
        inputs=[];spectra=[]
        for frame in range(3):
            t=np.arange(n)
            x=np.vstack([np.sin(2*np.pi*5*t/n)+.1*rng.normal(size=n)+.3,
                         .7*np.cos(2*np.pi*5*t/n)+.4*np.cos(2*np.pi*3*t/n)+.1*rng.normal(size=n)-.2])
            inputs.append(x)
            prepared=x if detrend=='None' else x-x.mean(axis=1,keepdims=True)
            estimates=[]
            for index in selected:
                taps=tapers[index]
                f=np.fft.rfft(prepared*taps,n=nfft,axis=1)/np.sqrt(128*np.sum(taps*taps))
                estimate=f.conj()[:,None,:]*f[None,:,:]
                estimate[:,:,1: -1 if nfft%2==0 else None]*=2
                estimates.append(estimate)
            spectra.append(np.mean(estimates,axis=0))
        result=np.mean(spectra,axis=0)
        cases.append(dict(length=n,nw=nw,candidates=k,floor=floor,nfft=nfft,rate=128.,detrend=detrend,
                          selected=selected.tolist(),concentrations=ratios.tolist(),tapers=tapers.tolist(),
                          inputs=np.asarray(inputs).tolist(),real=result.real.tolist(),imaginary=result.imag.tolist()))
    document=dict(schema='signal4s-multitaper-reference-v1',numpy=np.__version__,scipy=scipy.__version__,
                  tolerances=dict(taper_absolute=1e-9,concentration_absolute=1e-10,csd_relative=1e-9,csd_absolute=1e-12),
                  normalization='symmetric l2 DPSS; explicit concentration selection; equal taper mean then equal frame mean; one-sided density; conj(first)*second',cases=cases)
    payload=(json.dumps(document,indent=2)+'\n').encode()
    def vec(values):return 'Vector('+','.join(repr(float(v)) for v in values)+')'
    declarations=[]
    for i,c in enumerate(cases):
        fields=[str(c['length']),repr(c['nw']),str(c['candidates']),repr(c['floor']),str(c['nfft']),repr(c['rate']),json.dumps(c['detrend']),
              'Vector('+','.join(map(str,c['selected']))+')',vec(c['concentrations']),vec(np.array(c['tapers']).ravel()),
              vec(np.array(c['inputs']).ravel()),vec(np.array(c['real']).ravel()),vec(np.array(c['imaginary']).ravel())]
        declarations.append('  private def case'+str(i)+' = MultitaperReferenceCase('+','.join(fields)+')\n')
    scala='''package signal4s.fft
final case class MultitaperReferenceCase(length: Int,nw: Double,candidates: Int,floor: Double,nfft: Int,rate: Double,detrend: String,
  selected: Vector[Int],concentrations: Vector[Double],tapers: Vector[Double],inputs: Vector[Double],real: Vector[Double],imaginary: Vector[Double])
object MultitaperReference:
'''+''.join(declarations)+'  val cases = Vector('+','.join('case'+str(i) for i in range(len(cases)))+')\n'
    outputs={ROOT/'modules/fft/shared/src/test/resources/multitaper-reference-v1.json':payload,
             ROOT/'modules/fft/shared/src/test/scala/signal4s/fft/MultitaperReference.scala':scala.encode()}
    for path,data in outputs.items():
        if args.check:
            assert path.read_bytes()==data, str(path)
        else:
            assert not path.exists() or path.read_bytes()==data, 'refusing to replace different frozen reference'
            path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
        print(hashlib.sha256(data).hexdigest(),path.relative_to(ROOT))

if __name__=='__main__':main()
