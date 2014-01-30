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
      fact = dlog( x / xg(igrid,beta) ) / step(igrid)
      do j=0,beta-bound
         if(x.ge.xg(igrid,beta-j).and.x.lt.xg(igrid,beta-j+1))then
            w_int = 1d0
            do delta=0,k
               if(delta.ne.j) w_int = w_int * ( fact / ( j - delta ) 
     1                              + 1d0 )
            enddo
         endif
      enddo
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
      integer alpha,delta
      double precision step
      double precision xg(-10:n+10)
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
*     bound of the input grid and below the lower bound.
*
      step = dlog(xg(n)/xg(n-1))
      do alpha=n+1,n+k
         xg(alpha) = xg(alpha-1) * dexp(step)
      enddo
*
      step = dlog(xg(0)/xg(1))
      do alpha=-1,-k,-1
         xg(alpha) = xg(alpha+1) * dexp(step)
      enddo
*
      w_int_ext = 0d0
*     
      if(x.lt.xg(beta-k).or.x.ge.xg(beta+1)) return
*     
      do j=0,k
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
