package signal4s.multirate

import gale.linalg.{DVec, DVecBuilder}
import signal4s.SignalError

/** Type-I polyphase decomposition of a prototype FIR into `nPhases` branches.
  *
  * Phase `p` holds taps `h[p], h[p+nPhases], h[p+2*nPhases], …`.
  */
final case class PolyphaseBank private (
    phases: IArray[DVec],
    prototypeLength: Int
):
  def nPhases: Int = phases.length
  def maxPhaseLength: Int =
    var m = 0
    var i = 0
    while i < phases.length do
      m = math.max(m, phases(i).length)
      i += 1
    m

object PolyphaseBank:
  def decompose(h: DVec, nPhases: Int): Either[SignalError, PolyphaseBank] =
    if h.length == 0 then Left(SignalError.EmptyKernel)
    else if nPhases <= 0 then Left(SignalError.InvalidInputLength(nPhases))
    else
      val built = IArray.newBuilder[DVec]
      var p = 0
      while p < nPhases do
        val len = (h.length - 1 - p) / nPhases + (if p < h.length then 1 else 0)
        val taps = DVecBuilder.zeros(math.max(0, len))
        var i = 0
        var src = p
        while src < h.length do
          taps(i) = h(src)
          i += 1
          src += nPhases
        built += taps.result()
        p += 1
      Right(new PolyphaseBank(built.result(), h.length))
