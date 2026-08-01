# E6 FFT / OLA convolution bench

generated: 2026-08-01T22:57:22.009347Z

N=128 M=16 fftOneShotMs=0.1115 fftPlanReuseMs=0.0537 ola64OneShotMs=0.0584 ola64PlanReuseMs=0.0247
N=512 M=32 fftOneShotMs=0.0435 fftPlanReuseMs=0.0258 ola64OneShotMs=0.0324 ola64PlanReuseMs=0.0262
N=2048 M=64 fftOneShotMs=0.1728 fftPlanReuseMs=0.0881 ola64OneShotMs=0.0944 ola64PlanReuseMs=0.0948
N=8192 M=128 fftOneShotMs=0.3129 fftPlanReuseMs=0.1547 ola64OneShotMs=0.3178 ola64PlanReuseMs=0.2624

# Circular convolution

N=64 M=8 circularDirectMs=0.0076 circularFftMs=0.0050
N=128 M=16 circularDirectMs=0.0086 circularFftMs=0.0041
N=256 M=16 circularDirectMs=0.0079 circularFftMs=0.0064
N=512 M=32 circularDirectMs=0.0143 circularFftMs=0.0134
N=1024 M=48 circularDirectMs=0.0273 circularFftMs=0.0254
N=4096 M=80 circularDirectMs=0.0970 circularFftMs=0.1066
N=8192 M=128 circularDirectMs=0.2151 circularFftMs=0.1979

# Planned FFT region extraction

N=8192 M=128 region=valid fftPlanReuseMs=0.1306
N=8192 M=128 region=input fftPlanReuseMs=0.1323
N=16384 M=2048 region=valid fftPlanReuseMs=0.2878
N=16384 M=2048 region=input fftPlanReuseMs=0.2862
