************************************************************************
*
*     StructureFunctionxQ.f:
*
*     This function returns the value of the StructureFunctions at the
*     scale Q in GeV and for the bjorken variable x using the interpolation
*     of the cached StructureFunctions.
*
************************************************************************
      function StructureFunctionxQ(proc,sf,comp,x,Q)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/gridQ.h"
      include "../commons/scales.h"
      include "../commons/StructureFunctionsxQ.h"
      include "../commons/TimeLike.h"
**
*     Input Variables
*
      double precision x,Q
      character*2 proc
      character*2 sf
      character*6 comp
**
*     Internal Variables
*
      integer iproc,isf,icomp
      integer idx,idQ
      integer alpha,tau
      integer lx,ux,lQ,uQ
      integer isg,iQtrans,iQtransp
      integer tQ,ideg
      double precision Q2
      double precision w_int_xQ,wg
      double precision tol
      parameter(tol=1d-10)
**
*     Output Variables
*
      double precision StructureFunctionxQ
*
*     Check that the APFEL evolution has been initialized
*
      if(InCacheSFs.ne."done")then
         write(6,*) "In StructureFunctionxQ.f:"
         write(6,*) "Impossible to use this function because the"
         write(6,*) "structure functions have not been cached."
         write(6,*) "Call 'CacheStructureFunctionsAPFEL(Q0)'"
         write(6,*) "before calling StructureFunctionxQ."
         write(6,*) "   "
         call exit(-10)
      endif
*
      Q2 = Q * Q
*
*     Check consistency of the input variables
*
      if(proc.ne."EM".and.proc.ne."NC".and.proc.ne."CC")then
         write(6,*) "In StructureFunctionxQ.f:"
         write(6,*) "Invalid process, proc =",proc
         call exit(-10)
      endif
*
      if(proc.eq."EM") iproc = 0
      if(proc.eq."NC") iproc = 1
      if(proc.eq."CC") iproc = 2
*
      if(sf.ne."F2".and.sf.ne."FL".and.sf.ne."F3")then
         write(6,*) "In StructureFunctionxQ.f:"
         write(6,*) "Invalid StructureFunction, sf =",sf
         call exit(-10)
      endif
*
      if(sf.eq."F2") isf = 1
      if(sf.eq."FL") isf = 2
      if(sf.eq."F3") isf = 3
*
      if(comp(1:5).ne."light".and.comp(1:5).ne."charm".and.
     1   comp(1:6).ne."bottom".and.comp(1:3).ne."top".and.
     2   comp(1:5).ne."total")then
         write(6,*) "In StructureFunctionxQ.f:"
         write(6,*) "Invalid component, comp =",comp
         call exit(-10)
      endif
*
      if(comp(1:5).eq."light")  icomp = 3
      if(comp(1:5).eq."charm")  icomp = 4
      if(comp(1:6).eq."bottom") icomp = 5
      if(comp(1:3).eq."top")    icomp = 6
      if(comp(1:5).eq."total")  icomp = 7
*
*     Check that structure functions have been called in the allowed region
*
      if(x.lt.xmin(1)*(1d0-tol).or.x.gt.xmax*(1d0+tol))then
         write(6,*) "In StructureFunctionxQ.f:"
         write(6,*) "Value of x out of range, x =",x
         call exit(-10)
      endif
*
      if(Q2.lt.Q2min*(1d0-tol).or.Q2.gt.Q2max*(1d0+tol))then
         write(6,*) "In StructureFunctionxQ.f:"
         write(6,*) "Value of Q out of range, Q =",Q
         call exit(-10)
      endif
*
*     Trim away numerical fluctiations
*
      if(x.lt.xmin(1)) x  = xmin(1)
      if(x.gt.xmax)    x  = 1d0
*
      if(Q2.lt.Q2min)  Q2 = Q2min
      if(Q2.gt.Q2max)  Q2 = Q2max
*
*     1) If the value of Q2 happens to fall in the tiny gap
*     between two subsequent subgrids, push it on the first
*     node of the second grid.
*     2) Check that the Q2-subgrid has enough points for the
*     interpolation, otherwise reduce the interpolation degree.
*     3) Determine tQ, i.e. the parameter to avoid to interpolate
*     over the thresholds.
*
      idQ      = inter_degreeQ
      tQ       = 0
      iQtransp = 0
      iQtrans  = nQ(nfin) - 1
      do isg=nfin,nffi
*     1)
         if(Q2.gt.Q2g(iQtrans).and.Q2.lt.Q2g(iQtrans+1))
     1        Q2 = Q2g(iQtrans+1)
*     2)
         if(Q2.ge.Q2g(iQtransp).and.Q2.le.Q2g(iQtrans).and.
     1      nQ(isg).lt.idQ+1) idQ = nQ(isg) - 1
*     3)
         do ideg=2,idQ
            if(Q2.gt.Q2g(iQtrans-ideg+1).and.Q2.le.Q2g(iQtrans))
     1           tQ = ideg - 1
         enddo
         iQtransp = iQtrans + 1
         iQtrans  = iQtrans + nQ(isg+1)
      enddo
*
*     Interpolation
*
      idx = inter_degree(0)
      do alpha=0,nin(0)
         if(xg(0,alpha).gt.x) goto 11
      enddo
 11   lx = alpha - 1
      ux = lx + idx + 1
      do tau=0,nQ2g
         if(Q2g(tau).gt.Q2) goto 12
      enddo
 12   lQ = tau - 1 - tQ
      uQ = lQ + idQ + 1
      StructureFunctionxQ = 0d0
      do alpha=lx,ux
         do tau=lQ,uQ
            wg = w_int_xQ(tQ,idx,idQ,alpha,tau,x,Q2)
            StructureFunctionxQ = StructureFunctionxQ
     1           + wg * SFxQ(iproc,isf,icomp,alpha,tau)
         enddo
      enddo
*
      if(Timelike) StructureFunctionxQ = StructureFunctionxQ / x
*
      return
      end
