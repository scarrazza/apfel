************************************************************************
*
*     ExcludeScaleVariation.f:
*
*     This routine excludes from the precomputed coefficient functions
*     all the scale variation terms that were included in 
*     IncludeScaleVariation.f. This is needed when doing dynamical scale
*     variotions in such a way that any new iteration starts from a
*     "scale free" set of coefficient functions.
*
************************************************************************
      subroutine ExcludeScaleVariation
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
      include "../commons/MassInterpolIndices.h"
      include "../commons/mass_scheme.h"
      include "../commons/ColorFactors.h"
**
*     Internal Variables
*
      integer beta,alpha
      integer inf,jxi
      integer k
      integer gamma
      integer mapP(2,2),mapC(2)
      integer gbound
      integer ipt_FF
      double precision C12P0(4)
      double precision C1LP0(4)
      double precision C13P0(4)
      double precision P0P0(4)
      double precision tRen,tFac,tFac2h
      double precision beta0apf,b0,omlam
      double precision dh1,dlogxi
*
      if(ipt.eq.0) return
      if(krenQ.eq.1d0.and.kfacQ.eq.1d0) return
*
*     The minus signs serve to subtract the terms proportional to the logs
*
      tFac   = - dlog(kfacQ)
      tFac2h = - tFac * tFac / 2d0
      tRen   = tFac
*
*     Exclude renormalization scale variation terms due to the
*     running of the MSbar masses if needed
*
      dh1 = 0d0
      if(mass_scheme.eq."MSbar") dh1 = 6d0 * CF * tRen
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
*     NLO
*
         do beta=0,gbound
            do alpha=beta,nin(igrid)-1
*     Gluon
               SC2zm(igrid,inf,1,1,beta,alpha) = 
     1              SC2zm(igrid,inf,1,1,beta,alpha)
     2              - tFac * SP(igrid,inf,mapP(1,2),0,beta,alpha) / inf
               SC3zm(igrid,inf,1,1,beta,alpha) = 
     1              SC3zm(igrid,inf,1,1,beta,alpha)
     2              - tFac * SP(igrid,inf,mapP(1,2),0,beta,alpha) / inf
*     Plus
               SC2zm(igrid,inf,3,1,beta,alpha) = 
     1              SC2zm(igrid,inf,3,1,beta,alpha)
     2              - tFac * SP(igrid,inf,1,0,beta,alpha)
               SC3zm(igrid,inf,3,1,beta,alpha) = 
     1              SC3zm(igrid,inf,3,1,beta,alpha)
     2              - tFac * SP(igrid,inf,1,0,beta,alpha)
*     Minus
               SC2zm(igrid,inf,4,1,beta,alpha) = 
     1              SC2zm(igrid,inf,4,1,beta,alpha)
     2              - tFac * SP(igrid,inf,2,0,beta,alpha)
               SC3zm(igrid,inf,4,1,beta,alpha) = 
     1              SC3zm(igrid,inf,4,1,beta,alpha)
     2              - tFac * SP(igrid,inf,2,0,beta,alpha)
            enddo
         enddo
*
*     NNLO
*
         b0 = beta0apf(inf)
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
     2                 + tRen * b0 * SC2zm(igrid,inf,1,1,beta,alpha)
     3                 - tFac * C12P0(1)
     4                 + tFac2h * ( P0P0(1)
     5                 - b0 * SP(igrid,inf,mapP(1,2),0,beta,alpha) )
     6                 / inf
     7                 - tFac * SP(igrid,inf,mapP(1,2),1,beta,alpha)
     8                 / inf
                  SCLzm(igrid,inf,1,2,beta,alpha) = 
     1                 SCLzm(igrid,inf,1,2,beta,alpha)
     2                 + tRen * b0 * SCLzm(igrid,inf,1,1,beta,alpha)
     3                 - tFac * C1LP0(1)
                  SC3zm(igrid,inf,1,2,beta,alpha) = 
     1                 SC3zm(igrid,inf,1,2,beta,alpha)
     2                 + tRen * b0 * SC3zm(igrid,inf,1,1,beta,alpha)
     3                 - tFac * C13P0(1)
     4                 + tFac2h * ( P0P0(1)
     5                 - b0 * SP(igrid,inf,mapP(1,2),0,beta,alpha) )
     6                 / inf
     7                 - tFac * SP(igrid,inf,mapP(1,2),1,beta,alpha)
     8                 / inf
*     Pure-singlet
                  SC2zm(igrid,inf,2,2,beta,alpha) =
     1                 SC2zm(igrid,inf,2,2,beta,alpha)
     2                 - tFac * C12P0(2)
     3                 + tFac2h * P0P0(2) / inf
     4                 - tFac * ( SP(igrid,inf,mapP(1,1),1,beta,alpha)
     5                 - SP(igrid,inf,1,1,beta,alpha) ) / inf
                  SCLzm(igrid,inf,2,2,beta,alpha) =
     1                 SCLzm(igrid,inf,2,2,beta,alpha)
     2                 - tFac * C1LP0(2)
                  SC3zm(igrid,inf,2,2,beta,alpha) =
     1                 SC3zm(igrid,inf,2,2,beta,alpha)
     2                 - tFac * C13P0(2)
     3                 + tFac2h * P0P0(2) / inf
     3                 - tFac * ( SP(igrid,inf,mapP(1,1),1,beta,alpha)
     4                 - SP(igrid,inf,1,1,beta,alpha) ) / inf
*     Plus
                  SC2zm(igrid,inf,3,2,beta,alpha) = 
     1                 SC2zm(igrid,inf,3,2,beta,alpha)
     2                 + tRen * b0 * SC2zm(igrid,inf,3,1,beta,alpha)
     3                 - tFac * C12P0(3)
     4                 + tFac2h * ( P0P0(3)
     5                 - b0 * SP(igrid,inf,1,0,beta,alpha) )
     6                 - tFac * SP(igrid,inf,1,1,beta,alpha)
                  SCLzm(igrid,inf,3,2,beta,alpha) = 
     1                 SCLzm(igrid,inf,3,2,beta,alpha)
     2                 + tRen * b0 * SCLzm(igrid,inf,3,1,beta,alpha)
     3                 - tFac * C1LP0(3)
                  SC3zm(igrid,inf,3,2,beta,alpha) = 
     1                 SC3zm(igrid,inf,3,2,beta,alpha)
     2                 + tRen * b0 * SC3zm(igrid,inf,3,1,beta,alpha)
     3                 - tFac * C13P0(3)
     4                 + tFac2h * ( P0P0(3)
     5                 - b0 * SP(igrid,inf,1,0,beta,alpha) )
     6                 - tFac * SP(igrid,inf,1,1,beta,alpha)
*     Minus
                  SC2zm(igrid,inf,4,2,beta,alpha) = 
     1                 SC2zm(igrid,inf,4,2,beta,alpha)
     2                 + tRen * b0 * SC2zm(igrid,inf,4,1,beta,alpha)
     3                 - tFac * C12P0(4)
     4                 + tFac2h * ( P0P0(4)
     5                 - b0 * SP(igrid,inf,2,0,beta,alpha) )
     6                 - tFac * SP(igrid,inf,2,1,beta,alpha)
                  SCLzm(igrid,inf,4,2,beta,alpha) = 
     1                 SCLzm(igrid,inf,4,2,beta,alpha)
     2                 + tRen * b0 * SCLzm(igrid,inf,4,1,beta,alpha)
     3                 - tFac * C1LP0(4)
                  SC3zm(igrid,inf,4,2,beta,alpha) = 
     1                 SC3zm(igrid,inf,4,2,beta,alpha)
     2                 + tRen * b0 * SC3zm(igrid,inf,4,1,beta,alpha)
     3                 - tFac * C13P0(4)
     4                 + tFac2h * ( P0P0(4)
     5                 - b0 * SP(igrid,inf,2,0,beta,alpha) )
     6                 - tFac * SP(igrid,inf,2,1,beta,alpha)
               enddo
            enddo
         endif
      enddo
*
      ipt_FF = ipt
      if(MassScheme.eq."FONLL-B".and.ipt.ge.1) ipt_FF = 2
*
*     FFNS
*
      if(MassScheme(1:4).eq."FFNS".or.
     1   MassScheme(1:5).eq."FONLL")then
*
*     Neutral Current
*
         if(ipt_FF.ge.2)then
            b0 = beta0apf(Nf_FF)
            do beta=0,gbound
               do alpha=beta,nin(igrid)-1
                  do inf=4,6
                     do jxi=ixi(inf),ixi(inf)+1
                        if(ixi(inf).eq.0) cycle
*
*     Precompute needed convolutions
*
                        do k=1,2
                           C12P0(k) = 0d0
                           C1LP0(k) = 0d0
                        enddo
*
*     External grid
*
                        if(IsExt(igrid))then
                           do gamma=beta,alpha
*     Gluon
                              C12P0(1) = C12P0(1)
     4                         + SC2mNC(igrid,jxi,mapC(2),1,beta,gamma)
     5                         * SP(igrid,Nf_FF,mapP(2,2),0,gamma,alpha)
                              C1LP0(1) = C1LP0(1)
     4                         + SCLmNC(igrid,jxi,mapC(2),1,beta,gamma)
     5                         * SP(igrid,Nf_FF,mapP(2,2),0,gamma,alpha)
*     Pure-singlet
                              C12P0(2) = C12P0(2)
     1                         + SC2mNC(igrid,jxi,mapC(2),1,beta,gamma)
     2                         * SP(igrid,Nf_FF,mapP(2,1),0,gamma,alpha)
                              C1LP0(2) = C1LP0(2)
     1                         + SCLmNC(igrid,jxi,mapC(2),1,beta,gamma)
     2                         * SP(igrid,Nf_FF,mapP(2,1),0,gamma,alpha)
                           enddo
*
*     Internal grid
*
                        else
                           do gamma=beta,alpha
*     Gluon
                              C12P0(1) = C12P0(1)
     4                       + SC2mNC(igrid,jxi,mapC(2),1,0,gamma-beta)
     5                       * SP(igrid,Nf_FF,mapP(2,2),0,0,alpha-gamma)
                              C1LP0(1) = C1LP0(1)
     4                       + SCLmNC(igrid,jxi,mapC(2),1,0,gamma-beta)
     5                       * SP(igrid,Nf_FF,mapP(2,2),0,0,alpha-gamma)
*     Pure-singlet
                              C12P0(2) = C12P0(2)
     1                       + SC2mNC(igrid,jxi,mapC(2),1,0,gamma-beta)
     2                       * SP(igrid,Nf_FF,mapP(2,1),0,0,alpha-gamma)
                              C1LP0(2) = C1LP0(2)
     1                       + SCLmNC(igrid,jxi,mapC(2),1,0,gamma-beta)
     2                       * SP(igrid,Nf_FF,mapP(2,1),0,0,alpha-gamma)
                           enddo
                        endif
*
                        dlogxi = dlog( xigrid((jxi+1)*xistep)
     1                               / xigrid(jxi*xistep) )
                        if(jxi.eq.nxir) dlogxi = 1d8
*     Gluon
                        SC2mNC(igrid,jxi,1,2,beta,alpha) =
     1                       SC2mNC(igrid,jxi,1,2,beta,alpha)
     2                       + tRen * b0
     3                       * SC2mNC(igrid,jxi,1,1,beta,alpha)
     4                       - tFac * C12P0(1)
     5                       - dh1
     6                       * ( SC2mNC(igrid,jxi+1,1,1,beta,alpha) ! Numerical derivative
     7                       - SC2mNC(igrid,jxi,1,1,beta,alpha) )
     8                       / dlogxi
                        SCLmNC(igrid,jxi,1,2,beta,alpha) =
     1                       SCLmNC(igrid,jxi,1,2,beta,alpha)
     2                       + tRen * b0
     3                       * SCLmNC(igrid,jxi,1,1,beta,alpha)
     4                       - tFac * C1LP0(1)
     5                       - dh1
     6                       * ( SCLmNC(igrid,jxi+1,1,1,beta,alpha) ! Numerical derivative
     7                       - SCLmNC(igrid,jxi,1,1,beta,alpha) )
     8                       / dlogxi
*     Pure-singlet
                        SC2mNC(igrid,jxi,2,2,beta,alpha) =
     1                       SC2mNC(igrid,jxi,2,2,beta,alpha)
     2                       + tRen * b0
     3                       * SC2mNC(igrid,jxi,2,1,beta,alpha)
     4                       - tFac * C12P0(2)
                        SCLmNC(igrid,jxi,2,2,beta,alpha) =
     1                       SCLmNC(igrid,jxi,2,2,beta,alpha)
     2                       + tRen * b0
     3                       * SCLmNC(igrid,jxi,2,1,beta,alpha)
     4                       - tFac * C1LP0(2)
                     enddo
                  enddo
               enddo
            enddo
         endif
*
*     Charged Current
*     (Correct only NLO terms)
*
         do beta=0,gbound
            do alpha=beta,nin(igrid)-1
               do inf=4,6
                  do jxi=ixi(inf),ixi(inf)+1
                     if(ixi(inf).eq.0) cycle
                     omlam = 1d0 / ( 1d0 + xigrid(jxi*xistep) )
*     Gluon
                     SC2mCC(igrid,jxi,1,1,beta,alpha) = 
     1                    SC2mCC(igrid,jxi,1,1,beta,alpha)
     2                    - tFac
     3                    * SP(igrid,Nf_FF,mapP(1,2),0,beta,alpha)
     4                    / Nf_FF / 2d0
                     SCLmCC(igrid,jxi,1,1,beta,alpha) = 
     1                    SCLmCC(igrid,jxi,1,1,beta,alpha)
     2                    - tFac
     3                    * SP(igrid,Nf_FF,mapP(1,2),0,beta,alpha)
     4                    / Nf_FF  / 2d0 * omlam
                     SC3mCC(igrid,jxi,1,1,beta,alpha) = 
     1                    SC3mCC(igrid,jxi,1,1,beta,alpha)
     2                    - tFac
     3                    * SP(igrid,Nf_FF,mapP(1,2),0,beta,alpha)
     4                    / Nf_FF / 2d0
*     Plus
                     SC2mCC(igrid,jxi,3,1,beta,alpha) = 
     1                    SC2mCC(igrid,jxi,3,1,beta,alpha)
     2                    - tFac * SP(igrid,Nf_FF,1,0,beta,alpha)
                     SCLmCC(igrid,jxi,3,1,beta,alpha) = 
     1                    SCLmCC(igrid,jxi,3,1,beta,alpha)
     2                    - tFac * SP(igrid,Nf_FF,1,0,beta,alpha)
     3                    * omlam
                     SC3mCC(igrid,Nf_FF,3,1,beta,alpha) = 
     1                    SC3mCC(igrid,jxi,3,1,beta,alpha)
     2                    - tFac * SP(igrid,Nf_FF,1,0,beta,alpha)
                  enddo
               enddo
            enddo
         enddo
      endif
*
*     FFN0
*
      if(MassScheme(1:4).eq."FFN0".or.
     1   MassScheme(1:5).eq."FONLL")then
*
*     Neutral Current
*
         if(ipt_FF.ge.2)then
            b0 = beta0apf(Nf_FF)
            do beta=0,gbound
               do alpha=beta,nin(igrid)-1
                  do inf=4,6
                     do jxi=ixi(inf),ixi(inf)+1
                        if(ixi(inf).eq.0) cycle
*
*     Precompute needed convolutions
*
                        do k=1,2
                           C12P0(k) = 0d0
                           C1LP0(k) = 0d0
                        enddo
*
*     External grid
*
                        if(IsExt(igrid))then
                           do gamma=beta,alpha
*     Gluon
                              C12P0(1) = C12P0(1)
     1                         + SC2m0NC(igrid,jxi,mapC(2),1,beta,gamma)
     2                         * SP(igrid,Nf_FF,mapP(2,2),0,gamma,alpha)
                              C1LP0(1) = C1LP0(1)
     1                         + SCLm0NC(igrid,jxi,mapC(2),1,beta,gamma)
     2                         * SP(igrid,Nf_FF,mapP(2,2),0,gamma,alpha)
*     Pure-singlet
                              C12P0(2) = C12P0(2)
     1                         + SC2m0NC(igrid,jxi,mapC(2),1,beta,gamma)
     2                         * SP(igrid,Nf_FF,mapP(2,1),0,gamma,alpha)
                              C1LP0(2) = C1LP0(2)
     1                         + SCLm0NC(igrid,jxi,mapC(2),1,beta,gamma)
     2                         * SP(igrid,Nf_FF,mapP(2,1),0,gamma,alpha)
                           enddo
*
*     Internal grid
*
                        else
                           do gamma=beta,alpha
*     Gluon
                              C12P0(1) = C12P0(1)
     1                       + SC2m0NC(igrid,jxi,mapC(2),1,0,gamma-beta)
     2                       * SP(igrid,Nf_FF,mapP(2,2),0,0,alpha-gamma)
                              C1LP0(1) = C1LP0(1)
     1                       + SCLm0NC(igrid,jxi,mapC(2),1,0,gamma-beta)
     2                       * SP(igrid,Nf_FF,mapP(2,2),0,0,alpha-gamma)
*     Pure-singlet
                              C12P0(2) = C12P0(2)
     1                       + SC2m0NC(igrid,jxi,mapC(2),1,0,gamma-beta)
     2                       * SP(igrid,Nf_FF,mapP(2,1),0,0,alpha-gamma)
                              C1LP0(2) = C1LP0(2)
     1                       + SCLm0NC(igrid,jxi,mapC(2),1,0,gamma-beta)
     2                       * SP(igrid,Nf_FF,mapP(2,1),0,0,alpha-gamma)
                           enddo
                        endif
*
                        dlogxi = dlog( xigrid((jxi+1)*xistep)
     1                               / xigrid(jxi*xistep) )
                        if(jxi.eq.nxir) dlogxi = 1d8
*     Gluon
                        SC2m0NC(igrid,jxi,1,2,beta,alpha) =
     1                       SC2m0NC(igrid,jxi,1,2,beta,alpha)
     2                       + tRen * b0
     3                       * SC2m0NC(igrid,jxi,1,1,beta,alpha)
     4                       - tFac * C12P0(1)
     5                       - dh1
     6                       * ( SC2m0NC(igrid,jxi+1,1,1,beta,alpha) ! Numerical derivative
     7                       - SC2m0NC(igrid,jxi,1,1,beta,alpha) )
     8                       / dlogxi
                        SCLm0NC(igrid,jxi,1,2,beta,alpha) =
     1                       SCLm0NC(igrid,jxi,1,2,beta,alpha)
     2                       + tRen * b0
     3                       * SCLm0NC(igrid,jxi,1,1,beta,alpha)
     4                       - tFac * C1LP0(1)
*     Pure-singlet
                        SC2m0NC(igrid,jxi,2,2,beta,alpha) =
     1                       SC2m0NC(igrid,jxi,2,2,beta,alpha)
     2                       + tRen * b0
     3                       * SC2m0NC(igrid,jxi,2,1,beta,alpha)
     4                       - tFac * C12P0(2)
                        SCLm0NC(igrid,jxi,2,2,beta,alpha) =
     1                       SCLm0NC(igrid,jxi,2,2,beta,alpha)
     2                       + tRen * b0
     3                       * SCLm0NC(igrid,jxi,2,1,beta,alpha)
     4                       - tFac * C1LP0(2)
                     enddo
                  enddo
               enddo
            enddo
         endif
*
*     Charged Current
*     (Correct only NLO terms)
*
         do beta=0,gbound
            do alpha=beta,nin(igrid)-1
               do inf=4,6
                  do jxi=ixi(inf),ixi(inf)+1
                     if(ixi(inf).eq.0) cycle
*     Gluon
                     SC2m0CC(igrid,jxi,1,1,beta,alpha) = 
     1                    SC2m0CC(igrid,jxi,1,1,beta,alpha)
     2                    - tFac
     3                    * SP(igrid,Nf_FF,mapP(1,2),0,beta,alpha)
     4                    / Nf_FF / 2d0
                     SC3m0CC(igrid,jxi,1,1,beta,alpha) = 
     1                    SC3m0CC(igrid,jxi,1,1,beta,alpha)
     2                    - tFac
     3                    * SP(igrid,Nf_FF,mapP(1,2),0,beta,alpha)
     4                    / Nf_FF / 2d0
*     Plus
                     SC2m0CC(igrid,jxi,3,1,beta,alpha) = 
     1                    SC2m0CC(igrid,jxi,3,1,beta,alpha)
     2                    - tFac * SP(igrid,Nf_FF,1,0,beta,alpha)
                     SC3m0CC(igrid,Nf_FF,3,1,beta,alpha) = 
     1                    SC3m0CC(igrid,jxi,3,1,beta,alpha)
     2                    - tFac * SP(igrid,Nf_FF,1,0,beta,alpha)
                  enddo
               enddo
            enddo
         enddo
      endif
*
      return
      end
