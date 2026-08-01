package signal4s.laws

import signal4s.*
import signal4s.design.*

class DesignLawSuite extends munit.FunSuite:

  test("Butterworth ZPK/SOS/TF magnitude agreement"):
    val fs = SampleRate.hertz(1000.0).orThrow
    val d = Butterworth.lowPass(4, Frequency.hertz(100.0).orThrow, fs).orThrow
    DesignLaws.zpkSosTfAgree(d, fs)
    DesignLaws.butterworthStable(4, 100.0, 1000.0)

  test("highpass Butterworth stable"):
    DesignLaws.butterworthStable(3, 200.0, 1000.0)
    val fs = SampleRate.hertz(1000.0).orThrow
    val d = Butterworth.highPass(3, Frequency.hertz(200.0).orThrow, fs).orThrow
    DesignLaws.zpkSosTfAgree(d, fs, nfft = 32)
