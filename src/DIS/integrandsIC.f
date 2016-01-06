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
      double precision z,w_int,fR,fL
      double precision eta
      double precision CR,CL
      double precision c21ICR,cL1ICR,c31ICR
**
*     Output Variables
*
      double precision integrandsICm
*
      integrandsICm = 0d0
      if(y.ge.1d0) return
*
      eta = 2d0 * Q2IC / ( Spm + Del )
*
*     Interpolant functions
*
      z = xg(igrid,wbeta) / y / eta
*
      fR = w_int(inter_degree(igrid),walpha,z)
      fL = w_int(inter_degree(igrid),walpha,xg(igrid,wbeta)/eta)
*
*     Contructing integrands
*
*     C2
      if(sf.eq.1)then
         CR = c21ICR(y)
         CL = c21ICR(one)
*     CL
      elseif(sf.eq.2)then
         CR = cL1ICR(y)
         CL = cL1ICR(one)
*     C3
      elseif(sf.eq.3)then
         CR = c31ICR(y)
         CL = c31ICR(one)
      endif
*
      integrandsICm = ( CR * fR - CL * fL ) / ( 1d0 - y )
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
      double precision DICa,DICb
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
      fR = w_int(inter_degree(igrid),walpha,z)
      fS = fR - fL
*
      xi = xigrid(wixi)
*
*     Contructing integrands
*
      integrandsICm0 = DICa(xi,y) * fR + DICb(xi,y) * fS
*
      return
      end
