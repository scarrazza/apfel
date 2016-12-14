************************************************************************
*
*     xPDFxQ.f:
*
*     This function returns the value of the i-th PDF in the physical
*     basis at the scale Q in GeV and for the bjorken variable x using
*     the interpolation of the cached PDFs.
*
************************************************************************
      function xPDFxQ(id,x,Q)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/gridQ.h"
      include "../commons/scales.h"
      include "../commons/fphxQ.h"
**
*     Input Variables
*
      integer id
      double precision x,Q
**
*     Internal Variables
*
      integer i
      integer idx,idQ
      integer alpha,tau
      integer isg,iQtrans,iQtransp
      integer tQ,ideg
      double precision Q2
      double precision w_int_xQ,win(0:nint_max,0:nQ2g_max)
      double precision tol
      parameter(tol=1d-10)
**
*     Output Variables
*
      double precision xPDFxQ
*
*     Check that the APFEL evolution has been initialized
*
      if(InCachePDFs.ne."done")then
         write(6,*) "In xPDFxQ.f:"
         write(6,*) "Impossible to use this function because PDFs"
         write(6,*) "have not been cached."
         write(6,*) "Call 'CachePDFsAPFEL(Q0)' before calling xPDFxQ."
         write(6,*) "   "
         call exit(-10)
      endif
*
      Q2 = Q * Q
*
*     Check consistency of the input variables
*
      if(abs(id).le.6)then
         i = id
      elseif(id.eq.21)then
         i = 0
      elseif(abs(id).ne.11.or.abs(id).ne.13.or.abs(id).ne.15)then
         i = ( id - sign(1,id) * 9 ) / 2
      elseif(id.eq.22)then
         i = id
      else
         write(6,*) "In xPDFxQ.f:"
         write(6,*) "Invalid PDF index, id =",id
         call exit(-10)
      endif
*
*     Check that PDFs have been called in the allowed region
*
      if(x.lt.xmin(1)*(1d0-tol).or.x.gt.xmax*(1d0+tol))then
         write(6,*) "In xPDFxQ.f:"
         write(6,*) "Value of x out of range, x =",x
         call exit(-10)
      endif
*
      if(Q2.lt.Q2min*(1d0-tol).or.Q2.gt.Q2max*(1d0+tol))then
         write(6,*) "In xPDFxQ.f:"
         write(6,*) "Value of Q out of range, Q =",Q
         call exit(-10)
      endif
*
*     trim away numerical fluctiations
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
      xPDFxQ = 0d0
      idx = inter_degree(0)
*
*     Tabulate interpolation functions
*
      do tau=0,nQ2g
         do alpha=0,nin(0)
            win(alpha,tau) = w_int_xQ(tQ,idx,idQ,alpha,tau,x,Q2)
         enddo
      enddo
*
      if(abs(id).le.6)then
         do tau=0,nQ2g
            do alpha=0,nin(0)
               if(win(alpha,tau).eq.0d0) cycle
               xPDFxQ = xPDFxQ + win(alpha,tau) * fphxQ(i,alpha,tau)
            enddo
         enddo
      elseif(abs(id).eq.11.or.abs(id).eq.13.or.abs(id).eq.15)then
         do tau=0,nQ2g
            do alpha=0,nin(0)
               if(win(alpha,tau).eq.0d0) cycle
               xPDFxQ = xPDFxQ + win(alpha,tau) * fleptonxQ(i,alpha,tau)
            enddo
         enddo
      elseif(id.eq.22)then
         do tau=0,nQ2g
            do alpha=0,nin(0)
               if(win(alpha,tau).eq.0d0) cycle
               xPDFxQ = xPDFxQ + win(alpha,tau) * fgammaxQ(alpha,tau)
            enddo
         enddo
      endif
      if(dabs(xPDFxQ).le.1d-12) xPDFxQ = 0d0
*
      return
      end
*
************************************************************************
      subroutine xPDFxQall(x,Q,xf)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/gridQ.h"
      include "../commons/scales.h"
      include "../commons/fphxQ.h"
**
*     Input Variables
*
      double precision x,Q
**
*     Internal Variables
*
      integer i
      integer idx,idQ
      integer alpha,tau
      integer isg,iQtrans,iQtransp
      integer tQ,ideg
      double precision Q2
      double precision w_int_xQ,win(0:nint_max,0:nQ2g_max)
      double precision tol
      parameter(tol=1d-10)
**
*     Output Variables
*
      double precision xf(-6:6)
*
*     Check that the APFEL evolution has been initialized
*
      if(InCachePDFs.ne."done")then
         write(6,*) "In xPDFxQ.f:"
         write(6,*) "Impossible to use this function because PDFs"
         write(6,*) "have not been cached."
         write(6,*) "Call 'CachePDFsAPFEL(Q0)' before calling xPDFxQ."
         write(6,*) "   "
         call exit(-10)
      endif
*
      Q2 = Q * Q
*
*     Check that PDFs have been called in the allowed region
*
      if(x.lt.xmin(1)*(1d0-tol).or.x.gt.xmax*(1d0+tol))then
         write(6,*) "In xPDFxQ.f:"
         write(6,*) "Value of x out of range, x =",x
         call exit(-10)
      endif
*
      if(Q2.lt.Q2min*(1d0-tol).or.Q2.gt.Q2max*(1d0+tol))then
         write(6,*) "In xPDFxQ.f:"
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
      iQtrans  = nQ(nfin)
      do isg=nfin,nffi
*     1)
         if(Q2.gt.Q2g(iQtrans-1).and.Q2.lt.Q2g(iQtrans))then
            Q2 = Q2g(iQtrans)
         endif
*     2)
         if(Q2.ge.Q2g(iQtransp).and.Q2.lt.Q2g(iQtrans).and.
     1      nQ(nfin).lt.idQ+1) idQ = nQ(nfin) - 1
*     3)
         do ideg=2,idQ
            if(Q2.gt.Q2g(iQtrans-ideg).and.Q2.le.Q2g(iQtrans-1))
     1           tQ = ideg - 1
         enddo
         iQtransp = iQtransp + nQ(isg)
         iQtrans  = iQtrans  + nQ(isg+1)
      enddo
*
      idx = inter_degree(0)
*
*     Tabulate interpolation functions
*
      do tau=0,nQ2g
         do alpha=0,nin(0)
            win(alpha,tau) = w_int_xQ(tQ,idx,idQ,alpha,tau,x,Q2)
         enddo
      enddo
*
*     Interpolation
*
      do i=-6,6
         xf(i) = 0
         do tau=0,nQ2g
            do alpha=0,nin(0)
               if(win(alpha,tau).eq.0d0) cycle
               xf(i) = xf(i) + win(alpha,tau) * fphxQ(i,alpha,tau)
            enddo
         enddo
         if(dabs(xf(i)).le.1d-12) xf(i) = 0d0
      enddo
*
      return
      end
