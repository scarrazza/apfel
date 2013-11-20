************************************************************************
*
*     interpolants.f:
*
*     This routine return the alpha-th interpolant function of degree n
*     in the value x.
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
