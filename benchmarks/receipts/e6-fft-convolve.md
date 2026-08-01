# E6 FFT / OLA convolution bench

generated: 2026-08-01T22:09:49.374323Z

N=128 M=16 fftOneShotMs=0.1079 fftPlanReuseMs=0.0505 ola64OneShotMs=0.0675 ola64PlanReuseMs=0.0270
N=512 M=32 fftOneShotMs=0.0398 fftPlanReuseMs=0.0216 ola64OneShotMs=0.2746 ola64PlanReuseMs=0.1634
N=2048 M=64 fftOneShotMs=0.4543 fftPlanReuseMs=0.1573 ola64OneShotMs=0.1183 ola64PlanReuseMs=0.0917
N=8192 M=128 fftOneShotMs=0.3833 fftPlanReuseMs=0.1733 ola64OneShotMs=0.4019 ola64PlanReuseMs=0.3453

# Circular convolution

N=64 M=8 circularDirectMs=0.0396 circularFftMs=0.0065
N=128 M=16 circularDirectMs=0.0336 circularFftMs=0.0097
N=256 M=16 circularDirectMs=0.0350 circularFftMs=0.0080
N=512 M=32 circularDirectMs=0.1309 circularFftMs=0.0147
N=1024 M=48 circularDirectMs=0.5206 circularFftMs=0.0291
N=4096 M=80 circularDirectMs=0.3757 circularFftMs=0.2136
N=8192 M=128 circularDirectMs=0.9885 circularFftMs=0.1291
