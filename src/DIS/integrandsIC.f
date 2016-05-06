************************************************************************
*
*     integrandsIC.f:
*
*     This functions return the integrands need to compute the IC
*     contributions to the structure functions.
*     In order to wrap it in a form that can be fed to the integration 
*     function DGAUSS, it is an explicit function only of the integration 
*     variable y while it depends on the other variables:
*
*     1) the grid indices walpha and wbeta,
*     2) the structure function index:
*
*        sf  Structure Function
*     --------------------------
*        1          F2
*        2          FL
*        3          F3
*
*     that are contained in the common block wrapDIS.h
*
************************************************************************
*
*     IC massive coefficient functions
*
************************************************************************
      function integrandsICm(y)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/wrapIC.h"
      include "../commons/wrapDIS.h"
      include "../commons/coeffhqmellin.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR,fS,fL
      double precision xi
      double precision CR,CS
      double precision c21ICR,cL1ICR,c31ICR
      double precision DICa,DICb
      double precision X0QGA
**
*     Output Variables
*
      double precision integrandsICm
*
      integrandsICm = 0d0
      if(y.ge.one) return
*
*     Interpolant functions
*
      z = xg(igrid,wbeta) / y / eta
*
      fR = w_int(inter_degree(igrid),walpha,z)
      fL = w_int(inter_degree(igrid),walpha,xg(igrid,wbeta)/eta)
      fS = fR - fL
*
      xi  = xigrid(wixi)
*
*     Contructing integrands
*
      CR = 0d0
      CS = 0d0
*     C2
      if(sf.eq.1)then
*     Gluon
         if(k.eq.1)then
            CR = X0QGA(y,1)
*     Non-singlet
         elseif(k.eq.3.or.k.eq.4)then
            if(wl.eq.1)then
               CS = c21ICR(one) / ( 1d0 - y )
               CR = c21ICR(y) / ( 1d0 - y ) - CS
            elseif(wl.eq.2)then
               CR = DICa(xi,y)
               CS = DICb(xi,y)
            endif
         endif
*     CL
      elseif(sf.eq.2)then
*     Gluon
         if(k.eq.1)then
            CR = X0QGA(y,1)
*     Non-singlet
         elseif(k.eq.3.or.k.eq.4)then
            if(wl.eq.1)then
               CS = cL1ICR(one) / ( 1d0 - y )
               CR = cL1ICR(y) / ( 1d0 - y ) - CS
            elseif(wl.eq.2)then
               CR = DICa(xi,y)
               CS = DICb(xi,y)
            endif
         endif
*     C3
      elseif(sf.eq.3)then
*     Gluon
         if(k.eq.1)then
            CR = X0QGA(y,1)
*     Non-singlet
         elseif(k.eq.3.or.k.eq.4)then
            if(wl.eq.1)then
               CS = c31ICR(one) / ( 1d0 - y )
               CR = c31ICR(y) / ( 1d0 - y ) - CS
            elseif(wl.eq.2)then
               CR = DICa(xi,y)
               CS = DICb(xi,y)
            endif
         endif
      endif
*
      integrandsICm = CR * fR + CS * fS
*
      return
      end
*
************************************************************************
      function integrandsICm0(y)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/wrapDIS.h"
      include "../commons/coeffhqmellin.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR,fS,fL
      double precision xi
      double precision CR,CS
      double precision DICa,DICb
      double precision X0QGA
**
*     Output Variables
*
      double precision integrandsICm0
*
*     Interpolant functions
*
      z = xg(igrid,wbeta) / y
*
      fL = 0d0
      if(walpha.eq.wbeta) fL = 1d0
*
      xi = xigrid(wixi)
*
      fR = w_int(inter_degree(igrid),walpha,z)
      fS = fR - fL
*
*     Contructing integrands
*
      CR = 0d0
      CS = 0d0
*     Gluon
      if(k.eq.1)then
         CR = X0QGA(y,1)
*     Non-singlet
      elseif(k.eq.3.or.k.eq.4)then
         CR = DICa(xi,y)
         CS = DICb(xi,y)
      endif
*
      integrandsICm0 = CR * fR + CS * fS
*
      return
      end
