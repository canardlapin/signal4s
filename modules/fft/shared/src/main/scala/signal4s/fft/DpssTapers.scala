package signal4s.fft

import gale.linalg.{DVec,Matrix}
import gale.numeric.ExactSum
import gale.spectral.{Eigen,EigenSelection,EigenVectors,SpectralBackend}
import signal4s.SignalError

/** Symmetric, unit-energy DPSS ordered by the descending tridiagonal spectrum.
  * Only the requested leading candidates can be selected. Their actual band
  * concentrations, including rejected candidates, are retained for inspection.
  */
final class DpssTapers private (val length: Int,val timeBandwidth: Double,val candidateCount: Int,
    val minimumConcentration: Double,val candidateConcentrations: Vector[Double],
    val selectedIndices: Vector[Int],val taps: Vector[DVec],val constructionBytes: BigInt):
  val convention: WindowConvention=WindowConvention.Symmetric
  val retainedBytes: BigInt=BigInt(taps.size)*length*8+BigInt(candidateCount)*8+BigInt(taps.size)*4
  val count: Int=taps.size

object DpssTapers:
  val MaxLength: Int=512
  private[fft] def failure(detail: String): SignalError=SignalError.NumericalFailure("dpss",detail)

  /** Full portable symmetric eigendecomposition is deliberately bounded. */
  def inspect(length: Int,timeBandwidth: Double,candidates: Int,minimumConcentration: Double): Either[SignalError,BigInt] =
    if length<2 || length>MaxLength || !timeBandwidth.isFinite || timeBandwidth<=0 || timeBandwidth>=length/2.0 ||
      candidates<=0 || candidates>length || !minimumConcentration.isFinite || minimumConcentration<0 || minimumConcentration>1 then
      Left(failure("length 2..512, 0 < NW < length/2, 1..length candidates and concentration floor in [0,1] required"))
    else Right(BigInt(length)*length*96*8+BigInt(length)*candidates*16+
      BigInt(candidates)*16+ExactSum.resources.stateBytes+ExactSum.resources.ratioScratchBytes)

  def withinBudget(length: Int,timeBandwidth: Double,candidates: Int,minimumConcentration: Double,
      numericBudgetBytes: BigInt): Either[SignalError,DpssTapers] =
    for
      bytes <- inspect(length,timeBandwidth,candidates,minimumConcentration)
      _ <- if bytes<=numericBudgetBytes then Right(()) else Left(failure("construction exceeds declared primitive capacity"))
      result <- generate(length,timeBandwidth,candidates,minimumConcentration,bytes)
    yield result

  private def generate(n: Int,nw: Double,candidates: Int,floor: Double,bytes: BigInt): Either[SignalError,DpssTapers] =
    val bandwidth=nw/n
    val scale=n.toDouble*n
    // Slepian's commuting tridiagonal operator, divided by n^2 solely to
    // condition its eigensolve. The eigenvectors, and the bandwidth, are unchanged.
    val tridiagonal=Matrix.tabulate(n,n): (i,j) =>
      if i==j then math.pow((n-1-2.0*i)/2,2)*math.cos(2*math.Pi*bandwidth)/scale
      else if math.abs(i-j)==1 then
        val k=math.max(i,j).toDouble
        k*(n-k)/(2*scale)
      else 0.0
    given SpectralBackend=SpectralBackend.none
    Eigen.eigSymmetric(tridiagonal,EigenSelection.All,EigenVectors.Right).flatMap(_.requireConverged)
      .left.map(e=>failure(e.toString)).flatMap: spectrum =>
        val windows=Vector.tabulate(candidates): ordinal =>
          val column=n-1-ordinal
          // A deterministic sign convention; all auto/cross spectra are
          // invariant to a common taper sign. No l-infinity normalization.
          val first=(0 until n).find(i=>math.abs(spectrum.eigenvectors(i,column))>1e-12).getOrElse(0)
          val sign=if spectrum.eigenvectors(first,column)<0 then -1.0 else 1.0
          DVec.tabulate(n)(i=>sign*spectrum.eigenvectors(i,column))
        val orthogonality=(for i<-windows.indices;j<-0 to i yield
          val dot=(0 until n).foldLeft(0.0)((sum,k)=>sum+windows(i)(k)*windows(j)(k))
          math.abs(dot-(if i==j then 1.0 else 0.0))).max
        val ratios=Vector.newBuilder[Double]
        var ordinal=0
        var problem=if !orthogonality.isFinite || orthogonality>1e-9 then
          Some(failure("unit-energy/orthogonality certificate exceeds 1e-9")) else Option.empty[SignalError]
        while ordinal<candidates && problem.isEmpty do
          val v=windows(ordinal);val sum=ExactSum.zero()
          var i=0
          while i<n && problem.isEmpty do
            var j=0
            while j<n && problem.isEmpty do
              val lag=i-j
              val kernel=if lag==0 then 2*bandwidth else math.sin(2*math.Pi*bandwidth*lag)/(math.Pi*lag)
              sum.add((v(i)*v(j))*kernel) match
                case Left(e)=>problem=Some(failure(e.message))
                case Right(_)=>()
              j+=1
            i+=1
          val concentration=sum.value
          // The theoretical interval is [0,1]. Only this declared roundoff
          // envelope is clamped; larger excursions refuse the candidate set.
          val slack=64*math.ulp(1.0)*n
          if !concentration.isFinite || concentration < -slack || concentration>1+slack then
            problem=Some(failure("band concentration falls outside its numerical certificate"))
          else ratios+=concentration.max(0).min(1)
          ordinal+=1
        problem match
          case Some(error)=>Left(error)
          case None=>
            val concentrations=ratios.result()
            val kept=concentrations.indices.filter(i=>concentrations(i)>=floor).toVector
            if kept.isEmpty then Left(failure("no requested candidate meets the explicit concentration floor"))
            else Right(new DpssTapers(n,nw,candidates,floor,concentrations,kept,kept.map(windows),bytes))
