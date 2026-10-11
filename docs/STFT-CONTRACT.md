# STFT reconstruction coverage and frame coordinates

The default `StftPlan` prepares a canonical synthesis dual. Every hop residue
must have nonzero overlap energy. A periodic Hann window with hop equal to its
length has uncovered sample positions and is refused; hop greater than the
window length is also refused. An explicit custom synthesis window does not
establish this canonical invertibility guarantee.

Dual construction normalizes each hop residue separately, with no absolute
small-energy cutoff. Small but representable window scales remain supported.
Unrepresentable dual coefficients refuse. Numerical binary64 FFT/reconstruction
range limits still apply; this is not an exact-real arithmetic guarantee.

FramePlan's legacy Int count/start accessors use checked arithmetic and throw
ArithmeticException if their result cannot be represented. The additive
`frameCountChecked(Long)` and `frameStartChecked(Long)` accessors return typed
SignalError failures instead, use checked wide intermediates, and allocate no
input, output or FFT buffers. STFT analysis checks output frame-count capacity
before opening its transform workspace.

Time-axis labels refer to the discrete center sample:
`frameStart + floor(frameLength/2)`. A centered first frame therefore labels0;
a start-aligned first frame labels `floor(frameLength/2)/sampleRate` seconds.
Odd-length windows use the same discrete-center convention. Coefficient phase
remains the raw windowed rFFT convention; adjusting a label does not add a
centering phase to coefficients.

These corrections do not introduce streaming STFT or inverse guarantees for
the proposed bounded analysis-only API (GitHub issue19). The latter requires
separate workspace, normalization, ownership and consumer qualification.
