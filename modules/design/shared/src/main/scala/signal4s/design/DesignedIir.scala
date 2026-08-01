package signal4s.design

import signal4s.filter.SecondOrderCascade

/** Designed IIR with ZPK, SOS (primary for execution), and a report. */
final case class DesignedIir(
    zpk: ZerosPolesGain,
    sos: SecondOrderCascade,
    report: IirDesignReport
)
