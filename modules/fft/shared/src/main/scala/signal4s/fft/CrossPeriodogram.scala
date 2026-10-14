package signal4s.fft

import signal4s.SignalError

/** One shared convention for ordinary and multitaper common-window CSD. */
private[fft] object CrossPeriodogram:
  def write(spectra: Array[ComplexVector],nfft: Int,real: Array[Double],imaginary: Array[Double]): Either[SignalError,Unit] =
    val bins=nfft/2+1
    var index=0;var row=0
    while row<spectra.length do
      var column=row
      while column<spectra.length do
        var k=0
        while k<bins do
          val a=spectra(row);val b=spectra(column)
          val fold=if k==0 || (nfft%2==0 && k==bins-1) then 1.0 else 2.0
          val re=(a.real(k)*b.real(k)+a.imaginary(k)*b.imaginary(k))*fold
          val im=if row==column then 0.0 else (a.real(k)*b.imaginary(k)-a.imaginary(k)*b.real(k))*fold
          if !re.isFinite || !im.isFinite then return Left(BoundedCsdPlan.failure("nonfinite cross-periodogram"))
          real(index)=re;imaginary(index)=im
          index+=1;k+=1
        column+=1
      row+=1
    Right(())
