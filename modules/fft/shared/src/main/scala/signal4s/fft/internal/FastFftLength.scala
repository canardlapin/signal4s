package signal4s.fft.internal

/** Transform-length selection for real convolution FFTs.
  *
  * SciPy's `next_fast_len` can be shorter than the next power of two. Our portable
  * mixed-radix path for arbitrary 5-smooth lengths is correct but currently much
  * slower than the packed real power-of-two path, so [[apply]] prefers the next
  * power of two. The 5-smooth helpers remain for tests and for a future cost model
  * once mixed-radix is competitive—not for chasing a single benchmark cell.
  */
private[fft] object FastFftLength:

  def apply(minLength: Int): Int =
    require(minLength > 0, "minLength must be positive")
    nextPowerOfTwo(minLength)

  def isPowerOfTwo(n: Int): Boolean =
    n > 0 && (n & (n - 1)) == 0

  def isFiveSmooth(n: Int): Boolean =
    var x = n
    if x <= 0 then return false
    while x % 2 == 0 do x /= 2
    while x % 3 == 0 do x /= 3
    while x % 5 == 0 do x /= 5
    x == 1

  def nextPowerOfTwo(minLength: Int): Int =
    var n = 1
    while n < minLength do
      if n > Int.MaxValue / 2 then return Int.MaxValue
      n <<= 1
    n

  /** Smallest even 5-smooth integer `≥ minLength`. */
  def nextFiveSmooth(minLength: Int): Int =
    var n = minLength
    if n % 2 != 0 then n += 1
    while !isFiveSmooth(n) do
      if n >= Int.MaxValue - 1 then return Int.MaxValue
      n += 2
    n
