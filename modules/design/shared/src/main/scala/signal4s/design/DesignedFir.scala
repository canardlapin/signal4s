package signal4s.design

import signal4s.filter.Fir
import signal4s.fft.Window

final case class DesignedFir(
    fir: Fir,
    window: Window,
    band: FilterBand,
    numTaps: Int
)
