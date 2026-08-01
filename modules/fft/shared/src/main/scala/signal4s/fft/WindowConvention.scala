package signal4s.fft

/** DFT-periodicity convention for finite windows.
  *
  * Maps to SciPy `get_window(..., fftbins=)`:
  * [[Periodic]] ↔ `fftbins=True`, [[Symmetric]] ↔ `fftbins=False`.
  */
enum WindowConvention:
  case Symmetric
  case Periodic
