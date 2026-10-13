package signal4s.fft

import gale.linalg.{DVec, DVecBuilder}
import signal4s.{SampleRate, SignalError}

/** Discrete MNE-convention wavelet with dimensionless taps and norm sqrt(2). */
final class MorletKernel private[fft] (val frequency: Double, val cycles: Double,
    val halfWidth: Int, val taps: ComplexVector)

/** Owned frequency-by-time coefficients; invalid edge cells contain no estimate. */
final class MorletBlock private[fft] (val frequencies: Int, val samples: Int,
    val firstCenter: Int, val real: DVec, val imaginary: DVec,
    private val validity: Array[Boolean]):
  def isValid(frequency: Int, sample: Int): Boolean =
    require(frequency >= 0 && frequency < frequencies && sample >= 0 && sample < samples)
    validity(frequency * samples + sample)

/** Immutable finite-support kernels; no decimation or hidden edge padding. */
final class BoundedMorletPlan private (val sampleRate: SampleRate,
    val kernels: Vector[MorletKernel], val zeroMean: Boolean, val kernelBytes: BigInt):
  val maximumHalo: Int = kernels.map(_.halfWidth).max
  def peakBytes(samples: Int): BigInt = kernelBytes + BigInt(kernels.size) * samples * 33

  /** Analyze exactly these centers. Full kernel support is required for validity.
    * Caller input is read synchronously and never retained by the returned block.
    */
  def analyzeBlock(input: DVec, firstCenter: Int, samples: Int,
      numericBudgetBytes: BigInt): Either[SignalError, MorletBlock] =
    val cells = BigInt(kernels.size) * samples
    if samples < 0 || firstCenter < 0 || BigInt(firstCenter) + samples > input.length || cells > Int.MaxValue then
      return Left(BoundedMorletPlan.failure("centers or output shape outside input"))
    if peakBytes(samples) > numericBudgetBytes then return Left(BoundedMorletPlan.failure("coefficient output exceeds budget"))
    var i = 0
    while i < input.length do
      if !input(i).isFinite then return Left(BoundedMorletPlan.failure("nonfinite input"))
      i += 1
    val real = DVecBuilder.zeros(cells.toInt)
    val imaginary = DVecBuilder.zeros(cells.toInt)
    val valid = new Array[Boolean](cells.toInt)
    var frequency = 0
    while frequency < kernels.size do
      val kernel = kernels(frequency)
      var sample = 0
      while sample < samples do
        val center = firstCenter + sample
        if center >= kernel.halfWidth && center.toLong + kernel.halfWidth < input.length then
          var re = 0.0; var im = 0.0; var reCorrection = 0.0; var imCorrection = 0.0
          var tap = 0
          while tap < kernel.taps.length do
            val x = input(center + kernel.halfWidth - tap)
            val r = x * kernel.taps.real(tap) - reCorrection
            val nextRe = re + r
            reCorrection = (nextRe - re) - r
            re = nextRe
            val v = x * kernel.taps.imaginary(tap) - imCorrection
            val nextIm = im + v
            imCorrection = (nextIm - im) - v
            im = nextIm
            if !re.isFinite || !im.isFinite || !reCorrection.isFinite || !imCorrection.isFinite then
              return Left(BoundedMorletPlan.failure("convolution exceeds finite arithmetic capacity"))
            tap += 1
          val index = frequency * samples + sample
          real(index) = re; imaginary(index) = im; valid(index) = true
        sample += 1
      frequency += 1
    Right(new MorletBlock(kernels.size, samples, firstCenter, real.result(), imaginary.result(), valid))

object BoundedMorletPlan:
  val MaxKernelSamples: Int = 65535
  private[fft] def failure(message: String): SignalError = SignalError.NumericalFailure("morlet", message)

  def withinBudget(sampleRate: SampleRate, frequencies: Vector[Double], cycles: Vector[Double],
      numericBudgetBytes: BigInt, zeroMean: Boolean = true): Either[SignalError, BoundedMorletPlan] =
    val rate = sampleRate.hertz
    if !rate.isFinite || rate <= 0 || frequencies.isEmpty || !(cycles.size == 1 || cycles.size == frequencies.size) then
      return Left(failure("positive rate and nonempty matched frequency/cycle parameters required"))
    if !frequencies.forall(f => f.isFinite && f > 0 && f < rate / 2) ||
        !frequencies.sliding(2).forall(pair => pair.size < 2 || pair(0) < pair(1)) ||
        !cycles.forall(c => c.isFinite && c >= 1 && c <= 64) then
      return Left(failure("ordered sub-Nyquist frequencies and cycles in [1,64] required"))
    val dt = 1.0 / rate
    val sigmas = frequencies.zipWithIndex.map((f, i) => cycles(if cycles.size == 1 then 0 else i) / (2 * math.Pi * f))
    val halves = sigmas.map(sigma => math.ceil((5 * sigma) / dt) - 1)
    if !dt.isFinite || dt <= 0 || !sigmas.forall(s => s.isFinite && s > 0) ||
        !halves.forall(h => h.isFinite && h >= 0 && h <= (MaxKernelSamples - 1) / 2) then
      return Left(failure("finite support exceeds admitted kernel extent"))
    val bytes = halves.map(h => BigInt(h.toInt * 2 + 1) * 16).sum
    // Builders, detached taps and temporary unnormalized coefficients coexist.
    if bytes * 3 > numericBudgetBytes then return Left(failure("kernel construction exceeds budget"))
    val kernels = frequencies.indices.map: index =>
      val frequency = frequencies(index); val sigma = sigmas(index); val half = halves(index).toInt
      val length = half * 2 + 1
      val r = new Array[Double](length); val im = new Array[Double](length)
      val correction = if zeroMean then math.exp(-2 * math.pow(math.Pi * frequency * sigma, 2)) else 0.0
      var energy = 0.0; var compensation = 0.0; var j = 0
      while j < length do
        val time = (j - half) * dt
        val phase = 2 * math.Pi * frequency * time
        val envelope = math.exp(-(time * time) / (2 * sigma * sigma))
        r(j) = (math.cos(phase) - correction) * envelope
        im(j) = math.sin(phase) * envelope
        val term = r(j) * r(j) + im(j) * im(j) - compensation
        val next = energy + term
        compensation = (next - energy) - term
        energy = next
        j += 1
      val norm = math.sqrt(0.5) * math.sqrt(energy)
      new MorletKernel(frequency, cycles(if cycles.size == 1 then 0 else index), half,
        ComplexVector.unsafe(DVec.tabulate(length)(i => r(i) / norm), DVec.tabulate(length)(i => im(i) / norm)))
    // Refuse numerical construction failures before publishing any plan.
    if kernels.exists(k => !k.taps.real.toSeq.forall(_.isFinite) || !k.taps.imaginary.toSeq.forall(_.isFinite)) then
      Left(failure("degenerate or nonfinite wavelet normalization"))
    else Right(new BoundedMorletPlan(sampleRate, kernels.toVector, zeroMean, bytes))
