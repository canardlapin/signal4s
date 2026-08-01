package signal4s.fft.internal

/** Scala.js: no accelerated FFT backend. */
private[internal] object PlatformFft:
  def hasNativeRealBackend: Boolean = false
  def pow2(re: Array[Double], im: Array[Double], inverse: Boolean): Boolean = false
  def complex(re: Array[Double], im: Array[Double], inverse: Boolean): Boolean = false
  def supportsReal(n: Int): Boolean = false
  def realForward(re: Array[Double], im: Array[Double], onesided: Boolean): Boolean = false
  def realInverse(re: Array[Double], im: Array[Double]): Boolean = false
