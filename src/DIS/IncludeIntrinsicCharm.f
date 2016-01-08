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
**
*     Internal Variables
*
      integer ixi
      integer gbound
      integer alpha,beta
      integer bound
      double precision xi,lambda,eta
      double precision w_int,win
      double precision bq(0:6),dq(0:6),bqt(0:6)
      double precision fL
      double precision dgauss,c,d,eps
      double precision C2RS,C2L,CLRS,CLL,C3RS,C3L
      double precision c21ICL,cL1ICL,c31ICL
      double precision DIC,DICc
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
         xi     = xigrid(ixi*xistep)
         lambda = 1d0 / xi
         eta    = 2d0 / ( 1d0 + dsqrt( 1d0 + 4d0 * lambda ) )
*
*     LO
*
         do beta=0,gbound
            do alpha=beta,nin(igrid)-1
*
*     Neutral current
*
*     FFNS
               win = w_int(inter_degree(igrid),alpha,
     1              xg(igrid,beta)/eta)
*
               SC2mNC(igrid,ixi,3,0,beta,alpha) = ( 2d0 - eta ) * win
               SCLmNC(igrid,ixi,3,0,beta,alpha) = 4d0 * ( 1d0 - eta )
     1              / ( 2d0 - eta ) * win
*     FFN0
               if(alpha.eq.beta) SC2m0NC(igrid,ixi,3,0,beta,alpha) = 1d0
*
*     Charged current
*
*     (even if the IC contribution is a non-singlet one, put it in the 
*     pure-singlet slot (e.g. SC2mCC(igrid,ixi,2,0,beta,alpha) because
*     for now this slot is never used. If one day there will be the 
*     O(as^2) correction to CC, that contain a pure-singlet piece, I
*     will need to extend the CC arrays from 3 to 4.)
*
               if(alpha.eq.beta)then
*     FFNS
                  SC2mCC(igrid,ixi,2,0,beta,alpha) = ( 1d0 + lambda )
                  SCLmCC(igrid,ixi,2,0,beta,alpha) = lambda
                  SC3mCC(igrid,ixi,2,0,beta,alpha) = 1d0
*     FFN0
                  SC2m0CC(igrid,ixi,2,0,beta,alpha) = 1d0
                  SC3m0CC(igrid,ixi,2,0,beta,alpha) = 1d0
               endif
            enddo
         enddo
*
*     NLO
*
         if(ipt.ge.1)then
*
            wixi = ixi * xistep
*
*     FFNS
*
*     Charm mass and scale (Only pole mass for now)
            m12  = m2ph(4)
            Q2IC = m12 / lambda
*
*     Neutral current
*
*     Mass of the outcoming particle (charm)
            m22 = m12
*     Couplings
            call ComputeChargesDIS(Q2IC,bq,dq,bqt)
            Splus  = bq(4)
            Sminus = bqt(4)
            Rplus  = dq(4)        ! Needed only for F3 (never used)
            Rminus = 0d0          ! Needed only for F3 (never used)
*     Compute the needed IC factors according to hep-ph/9805233
            call ComputeICFactors("NC")
*
            do beta=0,gbound
               do alpha=beta,nin(igrid)-1
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
                  fL = w_int(inter_degree(igrid),alpha,
     1                 xg(igrid,beta)/eta)
*
                  walpha = alpha
                  wbeta  = beta
*
                  sf = 1
                  C2RS = eta * dgauss(integrandsICm,c,d,eps)
                  C2L  = eta * c21ICL(c)
*
                  sf = 2
                  CLRS = eta * dgauss(integrandsICm,c,d,eps)
                  CLL  = eta * cL1ICL(c)
*
                  SC2mNC(igrid,ixi,3,1,beta,alpha) = fact2
     1                                             * ( C2RS + C2L * fL )
                  SCLmNC(igrid,ixi,3,1,beta,alpha) = factL
     1                                             * ( CLRS + CLL * fL )
               enddo
            enddo
*
*     Charged current
*
*     Mass of the outcoming particle (strange or down (massless))
            m22 = m2strange
*     Couplings
            Splus  = 2d0 * ( V_cd2 + V_cs2 )
            Sminus = 0d0
            Rplus  = V_cd2 + V_cs2
            Rminus = 0d0
*     Compute the needed IC factors according to hep-ph/9805233
            call ComputeICFactors("CC")
*
            do beta=0,gbound
               do alpha=beta,nin(igrid)-1
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
                  fL = 0d0
                  if(alpha.eq.beta) fL = 1d0
*
                  walpha = alpha
                  wbeta  = beta
*
                  sf = 1
                  C2RS = dgauss(integrandsICm,c,d,eps)
                  C2L  = c21ICL(c)
*
                  sf = 2
                  CLRS = dgauss(integrandsICm,c,d,eps)
                  CLL  = cL1ICL(c)
*
                  sf = 3
                  C3RS = dgauss(integrandsICm,c,d,eps)
                  C3L  = c31ICL(c)
*
                  SC2mCC(igrid,ixi,2,1,beta,alpha) = fact2
     1                                             * ( C2RS + C2L * fL )
                  SCLmCC(igrid,ixi,2,1,beta,alpha) = factL
     1                                             * ( CLRS + CLL * fL )
                  SC3mCC(igrid,ixi,2,1,beta,alpha) = fact3
     1                                             * ( C3RS + C3L * fL )
               enddo
            enddo
*
*     FFN0
*
            do beta=0,gbound
               do alpha=beta,nin(igrid)-1
                  bound = alpha-inter_degree(igrid)
                  if(alpha.lt.inter_degree(igrid)) bound = 0
*
                  c = max(xg(igrid,beta),
     1                 xg(igrid,beta)/xg(igrid,alpha+1))
                  d = min(1d0,xg(igrid,beta)/xg(igrid,bound))
*
                  fL = 0d0
                  if(alpha.eq.beta) fL = 1d0
*
                  walpha = alpha
                  wbeta  = beta
*
                  DIC = dgauss(integrandsICm0,c,d,eps) + DICc(xi,c) * fL
*     Neutral Current
                  SC2m0NC(igrid,ixi,3,1,beta,alpha) =
     1                 SC2zm(igrid,Nf_FF,3,1,beta,alpha) + DIC
                  SCLm0NC(igrid,ixi,3,1,beta,alpha) =
     1                 SCLzm(igrid,Nf_FF,3,1,beta,alpha)
*     Charged Current
                  SC2m0CC(igrid,ixi,2,1,beta,alpha) =
     1                 SC2zm(igrid,Nf_FF,3,1,beta,alpha) + DIC
                  SCLm0CC(igrid,ixi,2,1,beta,alpha) =
     1                 SCLzm(igrid,Nf_FF,3,1,beta,alpha)
                  SC3m0CC(igrid,ixi,2,1,beta,alpha) =
     1                 SC3zm(igrid,Nf_FF,3,1,beta,alpha) + DIC
               enddo
            enddo
         endif
      enddo
*
      return
      end
*
************************************************************************
      subroutine ComputeICFactors(proc)
*
      implicit none
*
      include "../commons/wrapIC.h"
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      character*2 proc
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
*
      I1     = dlog( ( Spp + Del ) / ( Spp - Del ) ) / Del
      Cplus  = 2d0 * m1 * m2 * I1
      C1m    = - ( Spm * I1 + dlog( m12 / m22 ) ) / Q2IC
      C1p    = - ( Smp * I1 - dlog( m12 / m22 ) ) / Q2IC
      CRm    = ( Del2 / 2d0 / Q2IC
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
      if(proc.eq."NC")then
         factL = 2d0 * Splus / ( Splus + Sminus )
      elseif(proc.eq."CC")then
         factL = 1d0
      endif
*
      return
      end
