package signal4s.internal

/** Scala.js hot kernels for zero-boundary direct convolution. */
private[internal] object PlatformDirect:

  def preferFullForCausalInput(n: Int, m: Int): Boolean = false
  def preferFullForValid(n: Int, m: Int): Boolean = false

  def scatterFull(x: Array[Double], h: Array[Double], out: Array[Double]): Unit =
    val n = x.length
    val m = h.length
    val m4 = m & ~3
    var i = 0
    while i < n do
      val xi = x(i)
      if xi != 0.0 then
        var j = 0
        while j < m4 do
          out(i + j) += xi * h(j)
          out(i + j + 1) += xi * h(j + 1)
          out(i + j + 2) += xi * h(j + 2)
          out(i + j + 3) += xi * h(j + 3)
          j += 4
        while j < m do
          out(i + j) += xi * h(j)
          j += 1
      i += 1
