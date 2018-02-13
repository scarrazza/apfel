************************************************************************
*
*     IncludeIntrinsicCharm.f:
*
*     This routine includes in the precomputed coefficient functions
*     the intrinsic charm contributions to the massive coefficient
*     functions computed in RSLintegralsDIS.f.
*
************************************************************************
      subroutine IncludeIntrinsicCharm
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/m2th.h"
      include "../commons/grid.h"
      include "../commons/coeffhqmellin.h"
      include "../commons/integralsDIS.h"
      include "../commons/MassScheme.h"
      include "../commons/Nf_FF.h"
      include "../commons/CKM.h"
      include "../commons/wrapIC.h"
      include "../commons/wrapDIS.h"
      include "../commons/kfacQ.h"
**
*     Internal Variables
*
      integer ixi
      integer gbound
      integer alpha,beta
      integer bound
      double precision xi,lambda,lnF
      double precision w_int
      double precision bq(0:6),dq(0:6),bqt(0:6)
      double precision fL
      double precision dgauss,c,d,eps
      double precision C2RS,C2L,CLRS,CLL,C3RS,C3L
      double precision c21ICL,cL1ICL,c31ICL
      double precision DICc,gDIC,nsDIC
      double precision integrandsICm
      external integrandsICm
      double precision integrandsICm0
      external integrandsICm0
      parameter(eps=1d-5)
*
*     if the ZM-VFNS has been selected no IC contribution
*     has to be included.
*
      if(MassScheme.eq."ZM-VFNS") return
*
*     If an external grid is found compute the whole operator
*
      gbound = 0
      if(IsExt(igrid)) gbound = nin(igrid) - 1
*
      do ixi=1,nxir
         wixi   = ixi * xistep
*     Definitions
         xi     = xigrid(ixi*xistep)
         lambda = 1d0 / xi
         lnF    = dlog(xi*kfacQ)
*     Charm mass and scale.
*     Even though, due to the way how the coefficient functions are written
*     it is necessary to have mass and scale separarted, it can be shown
*     that the expressions only depend on the ratio.
         m12  = m2ph(4)
         Q2IC = m12 / lambda
*
*     FFNS
*
*     Neutral current
*
*     Mass of the outcoming particle (charm)
         m22 = m12
*     Couplings
         call ComputeChargesDIS(Q2IC,bq,dq,bqt)
         Splus  = bq(4)
         Sminus = bqt(4)
         Rplus  = dq(4)         ! Needed only for F3 (never used)
         Rminus = 0d0           ! Needed only for F3 (never used)
*     Compute the needed IC factors according to hep-ph/9805233
         call ComputeICFactors
*
         do beta=0,gbound
            do alpha=beta,nin(igrid)-1
*
               fL = w_int(inter_degree(igrid),alpha,
     1              xg(igrid,beta)/eta)
*
*     LO
*
               SC2mNC(igrid,ixi,3,0,beta,alpha) = ( 2d0 - eta ) * fL
               SCLmNC(igrid,ixi,3,0,beta,alpha) = 4d0 * ( 1d0 - eta )
     1              / ( 2d0 - eta ) / factL * fL
*
*     NLO
*
               if(ipt.ge.1)then
*
                  bound = alpha-inter_degree(igrid)
                  if(alpha.lt.inter_degree(igrid)) bound = 0
*
                  c = max(xg(igrid,beta),
     1                 xg(igrid,beta)/xg(igrid,alpha+1)) / eta
                  d = min(eta,xg(igrid,beta)/xg(igrid,bound)) / eta
*
                  if(c.ge.1d0) cycle
*
                  walpha = alpha
                  wbeta  = beta
*
                  k  = 3
                  wl = 1
*
                  sf   = 1
                  C2RS = fact2 * eta * dgauss(integrandsICm,c,d,eps)
                  C2L  = fact2 * eta * c21ICL(c)
*
                  sf   = 2
                  CLRS = eta * dgauss(integrandsICm,c,d,eps)
                  CLL  = eta * cL1ICL(c)
*
                  if(MassScheme(1:5).eq."FONLL")then
                     sf = 2
*     Gluon
                     k  = 1
                     gDIC = dgauss(integrandsICm,c,d,eps)
*     Put these terms in the LO gluon slot (e.g. SC2mNC(igrid,ixi,1,0,beta,alpha))
*     which is empty.
                     SC2mNC(igrid,ixi,1,0,beta,alpha) =
     1                    - ( 2d0 - eta ) * lnF * gDIC
                     SCLmNC(igrid,ixi,1,0,beta,alpha) =
     1                    - 4d0 * ( 1d0 - eta ) / ( 2d0 - eta ) / factL
     2                    * lnF * gDIC
*     Non-singlet
                     k  = 3
                     wl = 2
                     nsDIC = dgauss(integrandsICm,c,d,eps)
*
                     C2RS = C2RS - ( 2d0 - eta ) * nsDIC
                     C2L  = C2L - ( 2d0 - eta ) * DICc(xi,c)
*
                     CLRS = CLRS - 4d0 * ( 1d0 - eta ) / ( 2d0 - eta )
     1                    * nsDIC
                     CLL  = CLL - 4d0 * ( 1d0 - eta ) / ( 2d0 - eta )
     1                    * DICc(xi,c)
                  endif
*
                  SC2mNC(igrid,ixi,3,1,beta,alpha) = C2RS + C2L * fL
                  SCLmNC(igrid,ixi,3,1,beta,alpha) = CLRS + CLL * fL
               endif
            enddo
         enddo
*
*     Charged current
*
*     (even if the IC contribution is a non-singlet one, put it in the 
*     pure-singlet slot (e.g. SC2mCC(igrid,ixi,2,0,beta,alpha) because
*     for now this slot is never used. If one day there will be the 
*     O(as^2) correction to CC, that contain a pure-singlet piece, I
*     will need to extend the CC arrays from 3 to 4.)
*
*     Mass of the outcoming particle (strange or down (massless))
         m22 = m2strange
*     Couplings
         Splus  = 2d0 * ( V_cd2 + V_cs2 )
         Sminus = 0d0
         Rplus  = V_cd2 + V_cs2
         Rminus = 0d0
*     Compute the needed IC factors according to hep-ph/9805233
         call ComputeICFactors
*
         do beta=0,gbound
            do alpha=beta,nin(igrid)-1
*
               fL = 0d0
               if(alpha.eq.beta) fL = 1d0
*
*     LO
*
               SC2mCC(igrid,ixi,2,0,beta,alpha) = ( 1d0 + lambda ) * fL
               SCLmCC(igrid,ixi,2,0,beta,alpha) = lambda * fL
               SC3mCC(igrid,ixi,2,0,beta,alpha) = fL
*
*     NLO
*
               if(ipt.ge.1)then
*
                  bound = alpha-inter_degree(igrid)
                  if(alpha.lt.inter_degree(igrid)) bound = 0
*
                  c = max(xg(igrid,beta),
     1                    xg(igrid,beta)/xg(igrid,alpha+1))
                  d = min(1d0,xg(igrid,beta)/xg(igrid,bound))
*
                  if(c.ge.1d0) cycle
*
                  walpha = alpha
                  wbeta  = beta
*
                  k  = 3
                  wl = 1
*
                  sf   = 1
                  C2RS = fact2 * dgauss(integrandsICm,c,d,eps)
                  C2L  = fact2 * c21ICL(c)
*
                  sf   = 2
                  CLRS = dgauss(integrandsICm,c,d,eps)
                  CLL  = cL1ICL(c)
*
                  sf   = 3
                  C3RS = fact3 * dgauss(integrandsICm,c,d,eps)
                  C3L  = fact3 * c31ICL(c)
*
                  SC2mCC(igrid,ixi,2,1,beta,alpha) = C2RS + C2L * fL
                  SCLmCC(igrid,ixi,2,1,beta,alpha) = CLRS + CLL * fL
                  SC3mCC(igrid,ixi,2,1,beta,alpha) = C3RS + C3L * fL
               endif
            enddo
         enddo
*
*     FFN0
*
         do beta=0,gbound
            do alpha=beta,nin(igrid)-1
*
               fL = 0d0
               if(alpha.eq.beta) fL = 1d0
*
*     LO
*
*     Neutral Current
               SC2m0NC(igrid,ixi,3,0,beta,alpha) = fL
*     Charged Current
               SC2m0CC(igrid,ixi,2,0,beta,alpha) = fL
               SC3m0CC(igrid,ixi,2,0,beta,alpha) = fL
*
*     NLO
*
               if(ipt.ge.1)then
*
                  bound = alpha-inter_degree(igrid)
                  if(alpha.lt.inter_degree(igrid)) bound = 0
*
                  c = max(xg(igrid,beta),
     1                 xg(igrid,beta)/xg(igrid,alpha+1))
                  d = min(1d0,xg(igrid,beta)/xg(igrid,bound))
*
                  walpha = alpha
                  wbeta  = beta
*
                  k   = 3
                  nsDIC = dgauss(integrandsICm0,c,d,eps)
     1                 + DICc(xi,c) * fL
*     Neutral Current
                  SC2m0NC(igrid,ixi,3,1,beta,alpha) =
     1                 SC2zm(igrid,Nf_FF,3,1,beta,alpha) + nsDIC
                  SCLm0NC(igrid,ixi,3,1,beta,alpha) =
     1                 SCLzm(igrid,Nf_FF,3,1,beta,alpha)
*     Charged Current
                  SC2m0CC(igrid,ixi,2,1,beta,alpha) =
     1                 SC2zm(igrid,Nf_FF,3,1,beta,alpha) + nsDIC
                  SCLm0CC(igrid,ixi,2,1,beta,alpha) =
     1                 SCLzm(igrid,Nf_FF,3,1,beta,alpha)
                  SC3m0CC(igrid,ixi,2,1,beta,alpha) =
     1                 SC3zm(igrid,Nf_FF,3,1,beta,alpha) + nsDIC
*
*     If the FONLL scheme has been chosen, the non-singlet and the gluon
*     coefficient functions must be modified taking into account the
*     matching conditions. However, while the NC additional terms are
*     split into the FFNS and the FFN0 sector according to their kinematics,
*     the CC additional terms have all the same kinematics and thus I've
*     chosen to put them in the FFN0 terms.
*
                  if(MassScheme(1:5).eq."FONLL")then
*
*     Gluon
*
                     k = 1
                     gDIC = dgauss(integrandsICm0,c,d,eps)
*     Put these terms in the LO gluon slot which is empty.
*     Neutral Current
                     SC2m0NC(igrid,ixi,1,0,beta,alpha) = - lnF * gDIC
*     Charged Current
                     SC2m0CC(igrid,ixi,1,0,beta,alpha) =
     1                    lambda * lnF * gDIC / 2d0
                     SCLm0CC(igrid,ixi,1,0,beta,alpha) =
     1                    lambda * lnF * gDIC / 2d0
*
*     Non-singlet
*
*     Neutral Current
                     SC2m0NC(igrid,ixi,3,1,beta,alpha) =
     1                    SC2zm(igrid,Nf_FF,3,1,beta,alpha)
*     Charged Current
                     SC2m0CC(igrid,ixi,2,1,beta,alpha) =
     1                    SC2m0CC(igrid,ixi,2,1,beta,alpha)
     2                    + lambda * nsDIC
                     SCLm0CC(igrid,ixi,2,1,beta,alpha) =
     1                    SCLm0CC(igrid,ixi,2,1,beta,alpha)
     2                    + lambda * nsDIC
                  endif
               endif
            enddo
         enddo
      enddo
*
      return
      end
*
************************************************************************
      subroutine ComputeICFactors
*
      implicit none
*
      include "../commons/wrapIC.h"
      include "../commons/ColorFactors.h"
**
*     Internal Variables
*
      double precision DeltaFun
      double precision ddilog
*
      m1  = dsqrt(m12)
      m2  = dsqrt(m22)
*
      Del  = DeltaFun(m12,m22,-Q2IC)
      Del2 = Del * Del
      Spp  = Q2IC + m22 + m12
      Spm  = Q2IC + m22 - m12
      Smp  = Q2IC - m22 + m12
      eta  = 2d0 * Q2IC / ( Spm + Del )
*
      I1    = dlog( ( Spp + Del ) / ( Spp - Del ) ) / Del
      Cplus = 2d0 * m1 * m2 * I1
      C1m   = - ( Spm * I1 + dlog( m12 / m22 ) ) / Q2IC
      C1p   = - ( Smp * I1 - dlog( m12 / m22 ) ) / Q2IC
      CRm   = ( Del2 / 2d0 / Q2IC
     1     + Spp * ( 1d0 + dlog( Q2IC / Del ) ) ) * I1
     2     + ( m22 - m12 ) / 2d0 / Q2IC * dlog( m12 / m22 )
     3     - dlog( Q2IC / m12 ) - dlog( Q2IC / m22 ) - 4d0
     4     + Spp / Del * ( 
     5     + dlog( dabs( ( Del - Spm ) / 2d0 / Q2IC ) )**2 / 2d0
     6     + dlog( dabs( ( Del - Smp ) / 2d0 / Q2IC ) )**2 / 2d0
     7     - dlog( dabs( ( Del + Spm ) / 2d0 / Q2IC ) )**2 / 2d0
     8     - dlog( dabs( ( Del + Smp ) / 2d0 / Q2IC ) )**2 / 2d0
     9     - ddilog( ( Del - Spm ) / 2d0 / Del )
     1     - ddilog( ( Del - Smp ) / 2d0 / Del )
     2     + ddilog( ( Del + Spm ) / 2d0 / Del )
     3     + ddilog( ( Del + Smp ) / 2d0 / Del ) )
*
      S1 = 2d0 + Spp / Del * ( Del * I1
     1     + ddilog( 2d0 * Del / ( Del - Spp ) )
     2     - ddilog( 2d0 * Del / ( Del + Spp ) ) )
     3     + dlog( Del2 / m22 / Q2IC ) * ( - 2d0 + Spp * I1 )
      S2 = S1
      S3 = S1
*
      V1 = CRm
     1     + ( Sminus * Spp - 2d0 * Splus * m1 * m2 )
     2     / ( Splus * Spp - 2d0 * Sminus * m1 * m2 ) * Cplus
      V2 = CRm + ( m12 * C1p + m22 * C1m ) / 2d0 
     1     + Sminus / Splus * ( Cplus + m1 * m2 / 2d0
     2     * ( C1p + C1m ) )
      V3 = CRm + Rminus / Rplus * Cplus
*
      fact1 = 2d0 * CF * ( Spp - 2d0 * m1 * m2 * Sminus / Splus ) / Del
      fact2 = 2d0 * CF * Del / Q2IC
      fact3 = 2d0 * CF
      factL = 2d0 * Splus / ( Splus + Sminus )
*
      return
      end
