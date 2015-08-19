************************************************************************
*
*     IncludeScaleVariation.f:
*
*     This routine includes in the precomputed coefficient functions
*     all the scale variation terms that were not included in 
*     RSLintegralsDIS.f (or RSLintegralsDIS.f)
*
************************************************************************
      subroutine IncludeScaleVariation
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/grid.h"
      include "../commons/integrals.h"
      include "../commons/coeffhqmellin.h"
      include "../commons/integralsDIS.h"
      include "../commons/MassScheme.h"
      include "../commons/Nf_FF.h"
      include "../commons/krenQ.h"
      include "../commons/kfacQ.h"
**
*     Internal Variables
*
      integer beta,alpha
      integer inf,ixi
      integer k
      integer gamma
      integer mapP(2,2),mapC(2)
      integer gbound
      integer ipt_FF
      double precision C12P0(4)
      double precision C1LP0(4)
      double precision C13P0(4)
      double precision P0P0(4)
      double precision tR,tF,tF2h
      double precision beta0apf,b0
*
      if(ipt.eq.0) return
*
      tR   = dlog(krenQ)
      tF   = dlog(kfacQ)
      tF2h = tF * tF / 2d0
*
*     Maps used for the muliplications
*
      mapP(1,1) = 4
      mapP(1,2) = 5
      mapP(2,1) = 6
      mapP(2,2) = 7
*
      mapC(1) = 3
      mapC(2) = 1
*     
*     If an external grid is found compute the whole operator
*     
      gbound = 0
      if(IsExt(igrid)) gbound = nin(igrid) - 1
*
*     ZM-VFNS (Compute it always)
*
      do inf=3,6
*
         b0 = beta0apf(inf)
*
*     NNLO
*
         if(ipt.ge.2)then
            do beta=0,gbound
               do alpha=beta,nin(igrid)-1
*
*     Precompute needed convolutions
*
                  do k=1,4
                     P0P0(k)  = 0d0
                     C12P0(k) = 0d0
                     C1LP0(k) = 0d0
                     C13P0(k) = 0d0
                  enddo
*
*     External grid
*
                  if(IsExt(igrid))then
                     do gamma=beta,alpha
*     Gluon
                        P0P0(1) = P0P0(1)
     1                       + SP(igrid,inf,mapP(1,1),0,beta,gamma)
     2                       * SP(igrid,inf,mapP(1,2),0,gamma,alpha)
     3                       + SP(igrid,inf,mapP(1,2),0,beta,gamma)
     4                       * SP(igrid,inf,mapP(2,2),0,gamma,alpha)
                        C12P0(1) = C12P0(1)
     1                       + SC2zm(igrid,inf,mapC(1),1,beta,gamma)
     2                       * SP(igrid,inf,mapP(1,2),0,gamma,alpha)
     3                       / inf
     4                       + SC2zm(igrid,inf,mapC(2),1,beta,gamma)
     5                       * SP(igrid,inf,mapP(2,2),0,gamma,alpha)
                        C1LP0(1) = C1LP0(1)
     1                       + SCLzm(igrid,inf,mapC(1),1,beta,gamma)
     2                       * SP(igrid,inf,mapP(1,2),0,gamma,alpha)
     3                       / inf
     4                       + SCLzm(igrid,inf,mapC(2),1,beta,gamma)
     5                       * SP(igrid,inf,mapP(2,2),0,gamma,alpha)
                        C13P0(1) = C13P0(1)
     1                       + SC3zm(igrid,inf,mapC(1),1,beta,gamma)
     2                       * SP(igrid,inf,mapP(1,2),0,gamma,alpha)
     3                       / inf
     4                       + SC3zm(igrid,inf,mapC(2),1,beta,gamma)
     5                       * SP(igrid,inf,mapP(2,2),0,gamma,alpha)
*     Pure-singlet
                        P0P0(2) = P0P0(2)
     1                       + SP(igrid,inf,mapP(1,2),0,beta,gamma)
     2                       * SP(igrid,inf,mapP(2,1),0,gamma,alpha)
                        C12P0(2) = C12P0(2)
     1                       + SC2zm(igrid,inf,mapC(2),1,beta,gamma)
     2                       * SP(igrid,inf,mapP(2,1),0,gamma,alpha)
                        C1LP0(2) = C1LP0(2)
     1                       + SCLzm(igrid,inf,mapC(2),1,beta,gamma)
     2                       * SP(igrid,inf,mapP(2,1),0,gamma,alpha)
                        C13P0(2) = C13P0(2)
     1                       + SC3zm(igrid,inf,mapC(2),1,beta,gamma)
     2                       * SP(igrid,inf,mapP(2,1),0,gamma,alpha)
*     Plus
                        P0P0(3) = P0P0(3)
     1                       + SP(igrid,inf,1,0,beta,gamma)
     2                       * SP(igrid,inf,1,0,gamma,alpha)
                        C12P0(3) = C12P0(3)
     1                       + SC2zm(igrid,inf,3,1,beta,gamma)
     2                       * SP(igrid,inf,1,0,gamma,alpha)
                        C1LP0(3) = C1LP0(3)
     1                       + SCLzm(igrid,inf,3,1,beta,gamma)
     2                       * SP(igrid,inf,1,0,gamma,alpha)
                        C13P0(3) = C13P0(3)
     1                       + SC3zm(igrid,inf,3,1,beta,gamma)
     2                       * SP(igrid,inf,1,0,gamma,alpha)
                     enddo
*
*     Internal grid
*
                  else
                     do gamma=beta,alpha
*     Gluon
                        P0P0(1) = P0P0(1)
     1                       + SP(igrid,inf,mapP(1,1),0,0,gamma-beta)
     2                       * SP(igrid,inf,mapP(1,2),0,0,alpha-gamma)
     3                       + SP(igrid,inf,mapP(1,2),0,0,gamma-beta)
     4                       * SP(igrid,inf,mapP(2,2),0,0,alpha-gamma)
                        C12P0(1) = C12P0(1)
     1                       + SC2zm(igrid,inf,mapC(1),1,0,gamma-beta)
     2                       * SP(igrid,inf,mapP(1,2),0,0,alpha-gamma)
     3                       / inf
     4                       + SC2zm(igrid,inf,mapC(2),1,0,gamma-beta)
     5                       * SP(igrid,inf,mapP(2,2),0,0,alpha-gamma)
                        C1LP0(1) = C1LP0(1)
     1                       + SCLzm(igrid,inf,mapC(1),1,0,gamma-beta)
     2                       * SP(igrid,inf,mapP(1,2),0,0,alpha-gamma)
     3                       / inf
     4                       + SCLzm(igrid,inf,mapC(2),1,0,gamma-beta)
     5                       * SP(igrid,inf,mapP(2,2),0,0,alpha-gamma)
                        C13P0(1) = C13P0(1)
     1                       + SC3zm(igrid,inf,mapC(1),1,0,gamma-beta)
     2                       * SP(igrid,inf,mapP(1,2),0,0,alpha-gamma)
     3                       / inf
     4                       + SC3zm(igrid,inf,mapC(2),1,0,gamma-beta)
     5                       * SP(igrid,inf,mapP(2,2),0,0,alpha-gamma)
*     Pure-singlet
                        P0P0(2) = P0P0(2)
     1                       + SP(igrid,inf,mapP(1,2),0,0,gamma-beta)
     2                       * SP(igrid,inf,mapP(2,1),0,0,alpha-gamma)
                        C12P0(2) = C12P0(2)
     1                       + SC2zm(igrid,inf,mapC(2),1,0,gamma-beta)
     2                       * SP(igrid,inf,mapP(2,1),0,0,alpha-gamma)
                        C1LP0(2) = C1LP0(2)
     1                       + SCLzm(igrid,inf,mapC(2),1,0,gamma-beta)
     2                       * SP(igrid,inf,mapP(2,1),0,0,alpha-gamma)
                        C13P0(2) = C13P0(2)
     1                       + SC3zm(igrid,inf,mapC(2),1,0,gamma-beta)
     2                       * SP(igrid,inf,mapP(2,1),0,0,alpha-gamma)
*     Plus
                        P0P0(3) = P0P0(3)
     1                       + SP(igrid,inf,1,0,0,gamma-beta)
     2                       * SP(igrid,inf,1,0,0,alpha-gamma)
                        C12P0(3) = C12P0(3)
     1                       + SC2zm(igrid,inf,3,1,0,gamma-beta)
     2                       * SP(igrid,inf,1,0,0,alpha-gamma)
                        C1LP0(3) = C1LP0(3)
     1                       + SCLzm(igrid,inf,3,1,0,gamma-beta)
     2                       * SP(igrid,inf,1,0,0,alpha-gamma)
                        C13P0(3) = C13P0(3)
     1                       + SC3zm(igrid,inf,3,1,0,gamma-beta)
     2                       * SP(igrid,inf,1,0,0,alpha-gamma)
                     enddo
                  endif
*     Minus ( = Plus)
                  P0P0(4)  = P0P0(3)
                  C12P0(4) = C12P0(3)
                  C1LP0(4) = C1LP0(3)
                  C13P0(4) = C13P0(3)
*
*     Now adjust the NNLO coefficient functions
*
*     Gluon
                  SC2zm(igrid,inf,1,2,beta,alpha) = 
     1                 SC2zm(igrid,inf,1,2,beta,alpha)
     2                 + tR * b0 * SC2zm(igrid,inf,1,1,beta,alpha)
     3                 - tF * C12P0(1)
     4                 + tF2h * ( P0P0(1)
     5                 - b0 * SP(igrid,inf,mapP(1,2),0,beta,alpha) )
     6                 / inf
     7                 - tF * SP(igrid,inf,mapP(1,2),1,beta,alpha) / inf
                  SCLzm(igrid,inf,1,2,beta,alpha) = 
     1                 SCLzm(igrid,inf,1,2,beta,alpha)
     2                 + tR * b0 * SCLzm(igrid,inf,1,1,beta,alpha)
     3                 - tF * C1LP0(1)
                  SC3zm(igrid,inf,1,2,beta,alpha) = 
     1                 SC3zm(igrid,inf,1,2,beta,alpha)
     2                 + tR * b0 * SC3zm(igrid,inf,1,1,beta,alpha)
     3                 - tF * C13P0(1)
     4                 + tF2h * ( P0P0(1)
     5                 - b0 * SP(igrid,inf,mapP(1,2),0,beta,alpha) )
     6                 / inf
     7                 - tF * SP(igrid,inf,mapP(1,2),1,beta,alpha) / inf
*     Pure-singlet
                  SC2zm(igrid,inf,2,2,beta,alpha) =
     1                 SC2zm(igrid,inf,2,2,beta,alpha)
     2                 - tF * C12P0(2)
     3                 + tF2h * P0P0(2) / inf
     4                 - tF * ( SP(igrid,inf,mapP(1,1),1,beta,alpha)
     5                 - SP(igrid,inf,1,1,beta,alpha) ) / inf
                  SCLzm(igrid,inf,2,2,beta,alpha) =
     1                 SCLzm(igrid,inf,2,2,beta,alpha)
     2                 - tF * C1LP0(2)
                  SC3zm(igrid,inf,2,2,beta,alpha) =
     1                 SC3zm(igrid,inf,2,2,beta,alpha)
     2                 - tF * C13P0(2)
     3                 + tF2h * P0P0(2) / inf
     3                 - tF * ( SP(igrid,inf,mapP(1,1),1,beta,alpha)
     4                 - SP(igrid,inf,1,1,beta,alpha) ) / inf
*     Plus
                  SC2zm(igrid,inf,3,2,beta,alpha) = 
     1                 SC2zm(igrid,inf,3,2,beta,alpha)
     2                 + tR * b0 * SC2zm(igrid,inf,3,1,beta,alpha)
     3                 - tF * C12P0(3)
     4                 + tF2h * ( P0P0(3)
     5                 - b0 * SP(igrid,inf,1,0,beta,alpha) )
     6                 - tF * SP(igrid,inf,1,1,beta,alpha)
                  SCLzm(igrid,inf,3,2,beta,alpha) = 
     1                 SCLzm(igrid,inf,3,2,beta,alpha)
     2                 + tR * b0 * SCLzm(igrid,inf,3,1,beta,alpha)
     3                 - tF * C1LP0(3)
                  SC3zm(igrid,inf,3,2,beta,alpha) = 
     1                 SC3zm(igrid,inf,3,2,beta,alpha)
     2                 + tR * b0 * SC3zm(igrid,inf,3,1,beta,alpha)
     3                 - tF * C13P0(3)
     4                 + tF2h * ( P0P0(3)
     5                 - b0 * SP(igrid,inf,1,0,beta,alpha) )
     6                 - tF * SP(igrid,inf,1,1,beta,alpha)
*     Minus
                  SC2zm(igrid,inf,4,2,beta,alpha) = 
     1                 SC2zm(igrid,inf,4,2,beta,alpha)
     2                 + tR * b0 * SC2zm(igrid,inf,4,1,beta,alpha)
     3                 - tF * C12P0(4)
     4                 + tF2h * ( P0P0(4)
     5                 - b0 * SP(igrid,inf,2,0,beta,alpha) )
     6                 - tF * SP(igrid,inf,2,1,beta,alpha)
                  SCLzm(igrid,inf,4,2,beta,alpha) = 
     1                 SCLzm(igrid,inf,4,2,beta,alpha)
     2                 + tR * b0 * SCLzm(igrid,inf,4,1,beta,alpha)
     3                 - tF * C1LP0(4)
                  SC3zm(igrid,inf,4,2,beta,alpha) = 
     1                 SC3zm(igrid,inf,4,2,beta,alpha)
     2                 + tR * b0 * SC3zm(igrid,inf,4,1,beta,alpha)
     3                 - tF * C13P0(4)
     4                 + tF2h * ( P0P0(4)
     5                 - b0 * SP(igrid,inf,2,0,beta,alpha) )
     6                 - tF * SP(igrid,inf,2,1,beta,alpha)
               enddo
            enddo
         endif
*
*     NLO
*
         do beta=0,gbound
            do alpha=beta,nin(igrid)-1
*     Gluon
               SC2zm(igrid,inf,1,1,beta,alpha) = 
     1              SC2zm(igrid,inf,1,1,beta,alpha)
     2              - tF * SP(igrid,inf,mapP(1,2),0,beta,alpha) / inf
               SC3zm(igrid,inf,1,1,beta,alpha) = 
     1              SC3zm(igrid,inf,1,1,beta,alpha)
     2              - tF * SP(igrid,inf,mapP(1,2),0,beta,alpha) / inf
*     Plus
               SC2zm(igrid,inf,3,1,beta,alpha) = 
     1              SC2zm(igrid,inf,3,1,beta,alpha)
     2              - tF * SP(igrid,inf,1,0,beta,alpha)
               SC3zm(igrid,inf,3,1,beta,alpha) = 
     1              SC3zm(igrid,inf,3,1,beta,alpha)
     2              - tF * SP(igrid,inf,1,0,beta,alpha)
*     Minus
               SC2zm(igrid,inf,4,1,beta,alpha) = 
     1              SC2zm(igrid,inf,4,1,beta,alpha)
     2              - tF * SP(igrid,inf,2,0,beta,alpha)
               SC3zm(igrid,inf,4,1,beta,alpha) = 
     1              SC3zm(igrid,inf,4,1,beta,alpha)
     2              - tF * SP(igrid,inf,2,0,beta,alpha)
            enddo
         enddo
      enddo
*
      ipt_FF = ipt
      if(MassScheme.eq."FONLL-B".and.ipt.ge.1) ipt_FF = 2
*
*     FFNS
*
      if(MassScheme(1:4).eq."FFNS".or.
     1   MassScheme(1:5).eq."FONLL")then
         if(ipt_FF.ge.2)then
            b0 = beta0apf(Nf_FF)
            do beta=0,gbound
               do alpha=beta,nin(igrid)-1
                  do ixi=1,nxir
                     do k=1,3
                        SC2mNC(igrid,ixi,k,2,beta,alpha) = 
     1                       SC2mNC(igrid,ixi,k,2,beta,alpha)
     2                       + tR * b0
     3                       * SC2mNC(igrid,ixi,k,1,beta,alpha)
                        SCLmNC(igrid,ixi,k,2,beta,alpha) = 
     1                       SCLmNC(igrid,ixi,k,2,beta,alpha)
     2                       + tR * b0
     3                       * SCLmNC(igrid,ixi,k,1,beta,alpha)
                     enddo
                  enddo
               enddo
            enddo
         endif
      endif
*
*     FFN0
*
      if(MassScheme(1:4).eq."FFN0".or.
     1   MassScheme(1:5).eq."FONLL")then
         if(ipt_FF.ge.2)then
            b0 = beta0apf(Nf_FF)
            do beta=0,gbound
               do alpha=beta,nin(igrid)-1
                  do ixi=1,nxir
                     do k=1,3
                        SC2m0NC(igrid,ixi,k,2,beta,alpha) = 
     1                       SC2m0NC(igrid,ixi,k,2,beta,alpha)
     2                       + tR * b0
     3                       * SC2m0NC(igrid,ixi,k,1,beta,alpha)
                        SCLm0NC(igrid,ixi,k,2,beta,alpha) = 
     1                       SCLm0NC(igrid,ixi,k,2,beta,alpha)
     2                       + tR * b0
     3                       * SCLm0NC(igrid,ixi,k,1,beta,alpha)
                     enddo
                  enddo
               enddo
            enddo
         endif
      endif
*
      return
      end
