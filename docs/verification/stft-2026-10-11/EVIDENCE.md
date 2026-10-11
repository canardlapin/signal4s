# STFT coverage and frame-coordinate corrections

Baseline8821b1a, canonical main0b868a5 (same source tree). Linked ephys Motes ephys4s-provider-signal-stft-nola, ephys4s-provider-signal-frame-count and ephys4s-provider-signal-stft-time-axis; Signal4s GitHub20/21/22. Bounded analysis API19 remains separate.

Executed StftConsumerProbe reproduced a Hann4/hop4 default plan accepting a unit impulse at sample2 and reconstructing0; Rect4/hop5 incorrectly acquired a dual; FramePlan(4,2,Centered) returned-1073741824 instead of1073741824 for length Int.MaxValue-2; start-aligned Rect4/fs16 first center was labelled0 instead of0.125s. Baseline output and exact commands are preserved.

Canonical dual now requires nonzero coverage for every hop residue and refuses uncovered hop gaps. Per-residue normalization supports valid small/large/weak window scales and rejects unrepresentable dual coefficients. Frame count/start arithmetic is checked, with additive typed Long inspection methods and pre-workspace STFT output-count refusal. Labels use firstFrameStart+floor(windowLength/2), including odd windows, without altering coefficient phase.

Eight new regression cases pass on JVM and Scala.js. Targeted existing/new checks pass17/platform. Full FFT/STFT/Welch/reference checks pass141 JVM and107 JS, including the exact ephys4s SpectralContractsSuite12/platform injected through explicit unmanagedSources; STFT laws2/platform also pass. Zero failures/ignored checks. validation.json authenticates source, consumer, reference fixture and log hashes. Supported scopes remain separate from physiological accuracy or speed claims. Linux execution unavailable/skipped.

All provider JVM/JS builds were serialized with ephys4s comparator cohorts through its OS-held cohort lock. The measured ephys4s source code, classpath and provider lock were not changed; root pins remain old until the frozen schedule and its verification finish. This provider candidate uses a named independent branch, not an implicit consumer override.
