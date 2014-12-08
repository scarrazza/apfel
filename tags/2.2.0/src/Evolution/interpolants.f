************************************************************************
*
*     interpolants.f:
*
*     These routines return the beta-th interpolant function of degree
*     k in the value x.
*
*     "w_int" is used for the integrations. It is based on the assumption
*     of a logarithmic interpolation on logarithmically spaced grid but 
*     it's faster.
*
*     "w_int_gen" may be used for the final interpolation of PDFs to make
*     it continuous also at the junction between neighbour grids.
*     It is based on a logarithmic interpolation but on any grid.
*
*     "w_int_ext" same as "w_int_gen" but for a user-given grid different
*     from the internal grid of APFEL.
*
*     "w_int_herm" implements the Hermite cubic interpolation.
*
************************************************************************
      function w_int(k,beta,x)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer beta,k
      double precision x
**
*     Internal Variables
*
      integer j
      integer delta,bound
      double precision fact
**
*     Output Variables
*
      double precision w_int
*
      w_int = 0d0
      bound = beta - k
      if(k.gt.beta) bound = 0
      if(x.lt.xg(igrid,bound).or.x.ge.xg(igrid,beta+1)) return
*
      if(IsExt(igrid))then
         do j=0,beta-bound
            if(x.ge.xg(igrid,beta-j).and.x.lt.xg(igrid,beta-j+1))then
               w_int = 1d0
               do delta=0,k
                  if(delta.ne.j) w_int = w_int 
     1                 * dlog(x/xg(igrid,beta-j+delta)) 
     2                 / dlog(xg(igrid,beta)/xg(igrid,beta-j+delta))
               enddo
            endif
         enddo
      else
         fact = dlog( x / xg(igrid,beta) ) / step(igrid)
         do j=0,beta-bound
            if(x.ge.xg(igrid,beta-j).and.x.lt.xg(igrid,beta-j+1))then
               w_int = 1d0
               do delta=0,k
                  if(delta.ne.j) w_int = w_int * ( fact / ( j - delta ) 
     1                                 + 1d0 )
               enddo
            endif
         enddo
      endif
*
      return
      end
*
************************************************************************
      function w_int_gen(k,beta,x)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer beta,k
      double precision x
**
*     Internal Variables
*
      integer j
      integer delta,bound
**
*     Output Variables
*
      double precision w_int_gen
*
      w_int_gen = 0d0
      bound = beta - k
      if(k.gt.beta) bound = 0
      if(x.lt.xg(0,bound).or.x.ge.xg(0,beta+1)) return
*
      do j=0,beta-bound
         if(x.ge.xg(0,beta-j).and.x.lt.xg(0,beta-j+1))then
            w_int_gen = 1d0
            do delta=0,k
               if(delta.ne.j) w_int_gen = w_int_gen 
     1              * dlog(x/xg(0,beta-j+delta)) 
     2              / dlog(xg(0,beta)/xg(0,beta-j+delta))
            enddo
         endif
      enddo
*
      return
      end
*
************************************************************************
      function w_int_ext(n,xgi,k,beta,x)
*
      implicit none
**
*     Input Variables
*
      integer n,k,beta
      double precision xgi(0:n),x
**
*     Internal Variables
*
      integer j
      integer alpha,delta,bound
      double precision step
      double precision xg(0:n+10)
**
*     Output Variables
*
      double precision w_int_ext
*
*     Copy "xgi" into "xg"
*
      do alpha=0,n
         xg(alpha) = xgi(alpha)
      enddo
*
*     Create grid points beyond the upper
*     bound of the input grid.
*
      step = dlog(xg(n)/xg(n-1))
      do alpha=n+1,n+k
         xg(alpha) = xg(alpha-1) * dexp(step)
      enddo
*
      w_int_ext = 0d0
      bound = beta - k
      if(k.gt.beta) bound = 0
      if(x.lt.xg(bound).or.x.ge.xg(beta+1)) return
*
      do j=0,beta-bound
         if(x.ge.xg(beta-j).and.x.lt.xg(beta-j+1))then
            w_int_ext = 1d0
            do delta=0,k
               if(delta.ne.j) w_int_ext = w_int_ext 
     1              * dlog(x/xg(beta-j+delta)) 
     2              / dlog(xg(beta)/xg(beta-j+delta))
            enddo
         endif
      enddo
*     
      return
      end
*
************************************************************************
      function w_int_herm(n,xg,alpha,x)
*
      implicit none
**
*     Input Variables
*
      integer n,alpha
      double precision xg(0:200),x
**
*     Internal Variables
*
      integer k
      integer ub,lb
      double precision z,zg(-2:2),h(-2:1),t(-2:1)
      double precision h00,h10,h01,h11
      double precision Hp(-1:2)
      double precision eps
      parameter(eps=1d-14)
**
*     Output Variables
*
      double precision w_int_herm
*
*     Special case
*
      if(alpha.eq.n.and.dabs(x-1d0).lt.eps)then
         w_int_herm = 1d0
         return
      endif
*
      w_int_herm = 0d0
*
*     Return zero if out range
*
      lb = alpha - 2
      ub = alpha + 2
      if(lb.lt.0) lb = 0
      if(ub.gt.n) ub = n
      if(x.le.xg(lb)-eps.or.x.ge.xg(ub)+eps) return
*
*     Logarithmic interpolation
*
                         z      = dlog(x)
      if(alpha.ge.2)     zg(-2) = dlog(xg(alpha-2))
      if(alpha.ge.1)     zg(-1) = dlog(xg(alpha-1))
                         zg(0)  = dlog(xg(alpha))
      if(alpha.le.(n-1)) zg(1)  = dlog(xg(alpha+1))
      if(alpha.le.(n-2)) zg(2)  = dlog(xg(alpha+2))
c*
c*     Linear interpolation
c*
c                         z      = x
c      if(alpha.ge.2)     zg(-2) = xg(alpha-2)
c      if(alpha.ge.1)     zg(-1) = xg(alpha-1)
c                         zg(0)  = xg(alpha)
c      if(alpha.le.(n-1)) zg(1)  = xg(alpha+1)
c      if(alpha.le.(n-2)) zg(2)  = xg(alpha+2)
*
      h(-2) = zg(-1) - zg(-2)
      h(-1) = zg(0)  - zg(-1)
      h(0)  = zg(1)  - zg(0)
      h(1)  = zg(2)  - zg(1)
*
      t(-2) = ( z - zg(-2) ) / h(-2)
      t(-1) = ( z - zg(-1) ) / h(-1)
      t(0)  = ( z - zg(0) )  / h(0)
      t(1)  = ( z - zg(1) )  / h(1)
*     Hp(-1)
      if(alpha.ge.0.and.alpha.le.(n-2))then
         Hp(-1) = - h10(t(1)) * h(1) / h(0) / 2d0
      else
         Hp(-1) = 0d0
      endif
*     Hp(0)
      if(alpha.eq.0)then
         Hp(0)  = h00(t(0)) - h10(t(0)) - h11(t(0)) / 2d0
      elseif(alpha.ge.1.and.alpha.le.(n-2))then
         Hp(0)  = h00(t(0)) - h10(t(0)) * ( 1d0 - h(0) / h(-1) ) / 2d0
     1          - h11(t(0)) / 2d0
      elseif(alpha.eq.(n-1))then
         Hp(0)  = h00(t(0)) - h10(t(0)) * ( 1d0 - h(0) / h(-1) ) / 2d0
     1          - h11(t(0))
      elseif(alpha.eq.n)then
         Hp(0)  =   1
      endif
*     Hp(1)
      if(alpha.eq.0)then
         Hp(1)  = 0d0
      elseif(alpha.eq.1)then
         Hp(1)  = h01(t(-1)) + h11(t(-1)) * ( 1d0 - h(-1) / h(0) ) / 2d0
     1          + h10(t(-1))
      elseif(alpha.ge.2.and.alpha.le.(n-1))then
         Hp(1)  = h01(t(-1)) + h11(t(-1)) * ( 1d0 - h(-1) / h(0) ) / 2d0 
     1          + h10(t(-1)) / 2d0
      elseif(alpha.eq.n)then
         Hp(1)  = h01(t(-1)) + h11(t(-1)) + h10(t(-1)) / 2d0
      endif
*     Hp(2)
      if(alpha.ge.0.and.alpha.le.1)then
         Hp(2)  = 0d0
      else
         Hp(2)  = h11(t(-2)) * h(-2) / h(-1) / 2d0
      endif
*
      do k=-1,2
         if(x.ge.xg(alpha-k)-eps.and.x.lt.xg(alpha-k+1)-eps)
     1        w_int_herm = w_int_herm + Hp(k)
      enddo
*
      return
      end
*
************************************************************************
*
*     Hermite's polynomials
*
************************************************************************
      function h00(t)
*
      implicit none
**
*     Input Variables
*
      double precision t
**
*     Input Variables
*
      double precision h00
*
      h00 = 2d0 * t**3d0 - 3d0 * t**2d0 + 1d0
*
      return
      end
*
************************************************************************
      function h10(t)
*
      implicit none
**
*     Input Variables
*
      double precision t
**
*     Input Variables
*
      double precision h10
*
      h10 = t**3d0 - 2d0 * t**2d0 + t
*
      return
      end
*
************************************************************************
      function h01(t)
*
      implicit none
**
*     Input Variables
*
      double precision t
**
*     Input Variables
*
      double precision h01
*
      h01 = - 2d0 * t**3d0 + 3d0 * t**2d0
*
      return
      end
*
************************************************************************
      function h11(t)
*
      implicit none
**
*     Input Variables
*
      double precision t
**
*     Input Variables
*
      double precision h11
*
      h11 = t**3d0 - t**2d0
*
      return
      end
