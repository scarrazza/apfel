************************************************************************
*
*     sigmafk_dy.f:
*
*     Computates the FK tables.
*     It also computes the predictions unsing the internal APFEL PDFs
*     in two different ways:
*     - starting from the hard matrix coefficients at the final scale
*       using the APFEL evolution,
*     - using the FK tables at the initial scale using the PDFs at
*       the scale Q0.
*     The preditions should be in agreement.
*
************************************************************************
      subroutine sigmafk_dy(idat,Q0)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/mxdata.h"
      include "../commons/mxgridsize.h"
      include "../commons/mxgridsizeDY.h"
      include "../commons/xgridDY.h"
      include "../commons/DYcouplings.h"
      include "../commons/kinematics.h"
      include "../commons/cDY.h"
      include "../commons/sigmafk.h"
      include "../commons/xxDY.h"
      include "../commons/cutoff.h"
**
*     Input Variables
*
      integer idat
      double precision Q0
**
*     Internal Variables
*
      integer ixp(2)
      integer ibos
      integer ich
      integer nf
      integer ifl,ifl1,ifl2
      integer ipdf,jpdf
      integer ipdf1,ipdf2,jpdf1,jpdf2
      integer ix,jx
      integer ix1,ix2,jx1,jx2
      integer nIntervals,nx
      integer ipt,GetPerturbativeOrder
      integer i,kx
      integer ixpfk(2)
      double precision norm
      double precision zarat
      double precision tau,shad,m2,Q,as
      double precision gmu,GetGFermi
      double precision mz,GetZMass
      double precision mw,GetWMass
      double precision xGrid,xg(0:mxgridsize)
      double precision alphae
      double precision brZ
      double precision convfact
      double precision fact1,fact2
      double precision factorNS(-6:6,-6:6),factorQG(-6:6,-6:6)
      double precision sigmady1(-6:6,13,mxgridsizeDY,0:mxgridsize)
      double precision sigmady2(-6:6,13,mxgridsizeDY,0:mxgridsize)
      double precision lum_ns,lum_qg,lum_gq
      double precision ExternalEvolutionOperator
      double precision HeavyQuarkMass,a_QCD!,AlphaQCD
      double precision xPDF
      double precision f1Q(-6:6,mxgridsizeDY),f2Q(-6:6,mxgridsizeDY)
      double precision xfev(13),xfph(-6:6),xfQ0(13,0:mxgridsize)
      double precision pred_me,pred_fk,rel_diff
      double precision x0(2)
      double precision dist(0:mxgridsize)
      logical ComputePerdictions
      character*15 obslbl

      parameter(alphae = 7.297d-3)
      parameter(brZ = 0.033658d0)
      parameter(convfact = 0.389379304d9)
      parameter(ComputePerdictions = .true.)
*
*     Get obervable kinematics
*
      obslbl = obs(idat)
      m2     = q2dat(idat)
      tau    = x1dat(idat) * x2dat(idat)
      ixp(1) = ixp1DY(idat)
      ixp(2) = ixp2DY(idat)
      shad   = m2 / tau
      Q      = dsqrt(m2)
*
*     Physical constants
*
      gmu = GetGFermi()
      mz  = GetZMass()
      mw  = GetWMass()
*
*     Get perturbative order
*
      ipt = GetPerturbativeOrder()
*
*     Normalization and conversion factors
*
      if(obslbl.eq."DYP_E605")then
         zarat = 0.5d0      ! Isoscalar target
         ich   = 1
         ibos  = 1                    ! Photon production
         fact1 = ( 4d0 * pi * alphae**2d0 ) / ( 9d0 * m2 * shad )
         fact2 = 2d0 * shad**( 3d0 / 2d0 ) * dsqrt(m2)
         norm  = fact1 * fact2 * convfact * 1d-3
      elseif(obslbl.eq."DYP_E886P")then
         zarat = 1d0        ! Proton target
         ich   = 1
         ibos  = 1          ! Photon production
         fact1 = ( 4d0 * pi *alphae**2d0 ) / ( 9d0 * m2 * shad )
         fact2 = 2d0 * m2**2d0
         norm  = fact1 * fact2 * convfact * 1d-3
      elseif(obslbl.eq."DYP_E886D")then
         zarat = 0.5d0      ! Isoscalar target
         ich   = 1
         ibos  = 1          ! Photon production
         fact1 = ( 4d0 * pi *alphae**2d0 ) / ( 9d0 * m2 * shad )
         fact2 = 2d0 * m2**2d0
         norm  = fact1 * fact2 * convfact * 1d-3
      elseif(obslbl(1:9).eq."DYP_E906P")then
         zarat = 1d0        ! Proton target
         ich   = 1
         ibos  = 1          ! Photon production
         fact1 = ( 4d0 * pi *alphae**2d0 ) / ( 9d0 * m2 * shad )
         fact2 = 2d0 * m2**2d0
         norm  =  fact1 * fact2 * convfact * 1d-3
      elseif(obslbl(1:9).eq."DYP_E906D")then
         zarat = 0.5d0      ! Isoscalar target
         ich   = 1
         ibos  = 1              ! Photon production
         fact1 = ( 4d0 * pi *alphae**2d0 ) / ( 9d0 * m2 * shad )
         fact2 = 2d0 * m2**2d0
         norm  =  fact1 * fact2 * convfact * 1d-3
      elseif(obslbl.eq."EWK_ZRAP")then
         zarat = 1d0        ! antiproton target
         ich   = - 1
         ibos  = 2          ! Z production
         fact1 = pi * gmu * ( sqrt(2d0) / 3d0 ) * mz**2d0 / shad
         fact2 = 1d3 * brz
         norm  = fact1 * fact2 * convfact * 1d-3
      elseif(obslbl.eq."EWK_WASYM_WP")then
         zarat = 1d0        ! antiproton target
         ich   = - 1
         ibos  = 3          ! W+ production
         fact1 = pi * gmu * ( sqrt(2d0) / 3d0 ) * mw**2d0 / shad
         fact2 = 1d3
         norm  = fact1 * fact2 * convfact * 1d-3
      elseif(obslbl.eq."EWK_WASYM_WM")then
         zarat = 1d0        ! antiproton target
         ich   = - 1
         ibos  = 4          ! W- production
         fact1 = pi * gmu * ( sqrt(2d0) / 3d0 ) * mw**2d0 / shad
         fact2 = 1d3
         norm  = fact1 * fact2 * convfact * 1d-3
      elseif(obslbl(1:7).eq."DYP_PPY")then
         zarat = 1d0        ! proton target
         ich   = 1
         ibos  = 1          ! Photon production
         norm  = 1d0 !( 4d0 * pi * alphae**2d0 ) / ( 9d0 * m2 * shad )
      elseif(obslbl(1:7).eq."DYP_CPS")then
         zarat = 1d0        ! proton target
         ich   = 1
         ibos  = 3          ! W+ production
         norm  = 1d0
      elseif(obslbl(1:7).eq."DYP_CNS")then
         zarat = 1d0        ! proton target
         ich   = 1
         ibos  = 4          ! W- production
         norm  = 1d0
      elseif(obslbl(1:7).eq."DYP_CPO")then
         zarat = 1d0        ! antiproton target
         ich   = - 1
         ibos  = 3          ! W+ production
         norm  = 1d0
      elseif(obslbl(1:7).eq."DYP_CNO")then
         zarat = 1d0        ! antiproton target
         ich   = - 1
         ibos  = 4          ! W- production
         norm  = 1d0
      else
         write(6,*) "ERROR: in sigmafk_dy.f:"
         write(6,*) "Observable unsupported"
         write(6,*) "idat =",idat,"obslbl = ",obslbl 
         call exit(-10)         
      endif
*
*     Set number of active flavours
*
      if(Q.lt.HeavyQuarkMass(4,Q))then
         nf = 3
      elseif(Q.lt.HeavyQuarkMass(5,Q))then
         nf = 4
      elseif(Q.lt.HeavyQuarkMass(6,Q))then
         nf = 5
      else
         nf = 6
      endif
*
*     compute alphas
*
      as = 0d0
      !if(ipt.ge.1) as = AlphaQCD(Q) / 4d0 / pi
      if(ipt.ge.1) as = a_QCD(m2)
*
*     Initialize couplings
*
      do ipdf1=-nf,nf
         do ipdf2=-nf,nf
            factorNS(ipdf1,ipdf2) = 0d0
            factorQG(ipdf1,ipdf2) = 0d0
         enddo
      enddo
*
*     Calculate couplings. 
*     For the positivity constraint observables only the coupling of the 
*     correspondent flavour different from zero.
*
      if(obslbl(1:7).eq."DYP_PPY")then
         if(obslbl(13:14).eq."DW")then
            ifl = 1
         elseif(obslbl(13:14).eq."UP")then
            ifl = 2
         elseif(obslbl(13:14).eq."ST")then
            ifl = 3
         elseif(obslbl(13:14).eq."CH")then
            ifl = 4
         else
            write(6,*) "ERROR: in sigmafk_dy.f:"
            write(6,*) "Unknown flavour ",obslbl(13:14) 
            call exit(-10)
         endif
*
         factorNS(-ifl,-ifl) = CII(-ifl,ifl,ibos) 
     1        * ( VV(-ifl,ibos)**2d0 + AA(-ifl,ibos)**2d0 )
         factorNS(ifl,ifl)  = CII(ifl,-ifl,ibos) 
     1        * ( VV(ifl,ibos)**2d0 + AA(ifl,ibos)**2d0 )
*
         factorQG(-ifl,0) = CIF_NLO(nf,-ifl,ibos)
     1        * ( VV(-ifl,ibos)**2d0 + AA(-ifl,ibos)**2d0 )
         factorQG(0,-ifl) = factorQG(-ifl,0)
*
         factorQG(ifl,0)  = CIF_NLO(nf,ifl,ibos)
     1        * ( VV(ifl,ibos)**2d0 + AA(ifl,ibos)**2d0 )
         factorQG(0,ifl) = factorQG(ifl,0)
      elseif(obslbl(1:6).eq."DYP_CP")then
         if(obslbl(13:15).eq."UDB")then
            ifl1 = 2
            ifl2 = 1
         elseif(obslbl(13:15).eq."CDB")then
            ifl1 = 4
            ifl2 = 1
         elseif(obslbl(13:15).eq."TDB")then
            ifl1 = 6
            ifl2 = 1
         elseif(obslbl(13:15).eq."USB")then
            ifl1 = 2
            ifl2 = 3
         elseif(obslbl(13:15).eq."CSB")then
            ifl1 = 4
            ifl2 = 3
         elseif(obslbl(13:15).eq."TSB")then
            ifl1 = 6
            ifl2 = 3
         elseif(obslbl(13:15).eq."UBB")then
            ifl1 = 2
            ifl2 = 5
         elseif(obslbl(13:15).eq."CBB")then
            ifl1 = 4
            ifl2 = 5
         elseif(obslbl(13:15).eq."TBB")then
            ifl1 = 6
            ifl2 = 5
         else
            write(6,*) "ERROR: in sigmafk_dy.f:"
            write(6,*) "Unknown flavour ",obslbl(13:15) 
            call exit(-10)
         endif
*
*     Omit CKM matrix elements
*
         factorNS(ifl1,ifl2) = VV(ifl1,ibos)**2d0 + AA(ifl1,ibos)**2d0
         factorQG(ifl1,-ifl2) = factorNS(ifl1,ifl2)
      elseif(obslbl(1:6).eq."DYP_CN")then
         if(obslbl(13:15).eq."UBD")then
            ifl1 = - 2
            ifl2 = - 1
         elseif(obslbl(13:15).eq."CBD")then
            ifl1 = - 4
            ifl2 = - 1
         elseif(obslbl(13:15).eq."TBD")then
            ifl1 = - 6
            ifl2 = - 1
         elseif(obslbl(13:15).eq."UBS")then
            ifl1 = - 2
            ifl2 = - 3
         elseif(obslbl(13:15).eq."CBS")then
            ifl1 = - 4
            ifl2 = - 3
         elseif(obslbl(13:15).eq."TBS")then
            ifl1 = - 6
            ifl2 = - 3
         elseif(obslbl(13:15).eq."UBB")then
            ifl1 = - 2
            ifl2 = - 5
         elseif(obslbl(13:15).eq."CBB")then
            ifl1 = - 4
            ifl2 = - 5
         elseif(obslbl(13:15).eq."TBB")then
            ifl1 = - 6
            ifl2 = - 5
         else
            write(6,*) "ERROR: in sigmafk_dy.f:"
            write(6,*) "Unknown flavour ",obslbl(13:15) 
            call exit(-10)
         endif
*
*     Omit CKM matrix elements
*
         factorNS(ifl1,ifl2) = VV(ifl1,ibos)**2d0 + AA(ifl1,ibos)**2d0
         factorQG(ifl1,-ifl2) = factorNS(ifl1,ifl2)
      elseif(obslbl.eq."EWK_WASYM")then
         do ipdf1=-nf,nf
            do ipdf2=-nf,nf
               factorNS(ipdf1,ipdf2) = CII(ipdf1,-ipdf2,ibos)
     1         * ( VV(ipdf1,ibos)**2d0 + AA(ipdf1,ibos)**2d0 )
               factorQG(ipdf1,ipdf2) = 
     1         ( CIF(ipdf1,ipdf2,3) + CIF(ipdf1,ipdf2,4) )
     2         * ( VV(ipdf1,ibos)**2d0 + AA(ipdf1,ibos)**2d0 )/2d0
            enddo
         enddo
      else
         do ipdf=-nf,nf
           factorNS(ipdf,ipdf) = CII(ipdf,-ipdf,ibos)
     1     * ( VV(ipdf,ibos)**2d0 + AA(ipdf,ibos)**2d0 )
           factorQG(ipdf,0) = CIF_NLO(nf,ipdf,ibos)
     1     * ( VV(ipdf,ibos)**2d0 + AA(ipdf,ibos)**2d0 )
           factorQG(0,ipdf) = factorQG(ipdf,0)
         enddo
      endif
*
*     Perform PDF evolution with APFEL
*
      call EvolveAPFEL(Q0,Q)
*
*     Evaluate evolution matices
*
      nx = nIntervals()
      do jx=0,nx
         xg(jx) = xGrid(jx)
      enddo
*
      do ix=1,nxDY
         do jx=0,nx
            if(xxDY(ix).gt.xg(min(jx+1,nx))) cycle
            do jpdf=1,13
*     Projectile evolution matrix (always a proton)
               do ipdf=-6,6
                  sigmady1(ipdf,jpdf,ix,jx) = 
     1          ExternalEvolutionOperator("Ev2Ph",ipdf,jpdf,xxDY(ix),jx)
     2          / xxDY(ix)
               enddo
*     Target evolution matrix
               sigmady2(0,jpdf,ix,jx)  = sigmady1(0,jpdf,ix,jx)
               sigmady2(1,jpdf,ix,jx)  = 
     1              zarat * sigmady1(ich,jpdf,ix,jx)
     2              + ( 1d0 - zarat ) * sigmady1(ich*2,jpdf,ix,jx)
               sigmady2(2,jpdf,ix,jx)  = 
     1              ( 1d0 - zarat ) * sigmady1(ich,jpdf,ix,jx)
     2              + zarat * sigmady1(ich*2,jpdf,ix,jx)
               sigmady2(-1,jpdf,ix,jx) = 
     1              zarat * sigmady1(-ich,jpdf,ix,jx)
     2              + ( 1d0 - zarat ) * sigmady1(-ich*2,jpdf,ix,jx)
               sigmady2(-2,jpdf,ix,jx) =
     1              ( 1d0 - zarat ) * sigmady1(-ich,jpdf,ix,jx)
     2              + zarat * sigmady1(-ich*2,jpdf,ix,jx)
               do ipdf=3,6
                  sigmady2(ipdf,jpdf,ix,jx)  = 
     1                 sigmady1(ich*ipdf,jpdf,ix,jx)
                  sigmady2(-ipdf,jpdf,ix,jx) = 
     1                 sigmady1(-ich*ipdf,jpdf,ix,jx)
               enddo
            enddo
         enddo
      enddo
*
*     Generating FK grid
*
*     Initialization
*
      do jpdf1=1,13
         do jpdf2=1,13
            do jx1=0,nx
               do jx2=0,nx
                  sigmafkdy(jx1,jx2,jpdf1,jpdf2) = 0d0
               enddo
            enddo
         enddo
      enddo
*     Loop on the evolution basis PDFs
      do jpdf1=1,13
         do jpdf2=1,13
*     Loop on the Hard cross section grid points
            do ix1=ixp(1),nxDY
               do ix2=ixp(2),nxDY
*     Loop on the PDF luminosity grid points
                  do jx1=0,nx
                     if(xxDY(ix1).gt.xg(min(jx1+1,nx))) cycle
                     do jx2=0,nx
                        if(xxDY(ix2).gt.xg(min(jx2+1,nx))) cycle
*     Initialize luminosities
                        lum_ns = 0d0
                        lum_qg = 0d0
                        lum_gq = 0d0
*     W^+/W^- observables
                        if(obslbl(1:9).eq."EWK_WASYM".or.
     1                     obslbl(1:5).eq."DYP_C")then
                           do ipdf1=-nf,nf
                              do ipdf2=-nf,nf
*     quark-antiquark luminosity
                                 lum_ns = lum_ns
     1                                  + factorNS(ipdf1,ipdf2)
     2                                  * sigmady1(ipdf1,jpdf1,ix1,jx1)
     3                                  * sigmady2(-ipdf2,jpdf2,ix2,jx2)
*     quark-gluon luminosity
                                 lum_qg = lum_qg
     1                                  + factorQG(ipdf1,ipdf2)
     2                                  * sigmady1(ipdf1,jpdf1,ix1,jx1)
     3                                  * sigmady2(0,jpdf2,ix2,jx2)    
*     gluon-quark luminosity
                                 lum_gq = lum_gq
     1                                  + factorQG(ipdf1,ipdf2)
     2                                  * sigmady1(0,jpdf1,ix1,jx1)
     3                                  * sigmady2(ipdf2,jpdf2,ix2,jx2)
                              enddo
                           enddo
*     gamma/Z observables
                        else
                           do ipdf1=-nf,nf
                              lum_ns = lum_ns
     1                               + factorNS(ipdf1,ipdf1)
     2                               * sigmady1(ipdf1,jpdf1,ix1,jx1)
     3                               * sigmady2(-ipdf1,jpdf2,ix2,jx2)
                              lum_qg = lum_qg
     1                               + factorQG(ipdf1,0)
     2                               * sigmady1(ipdf1,jpdf1,ix1,jx1)
     3                               * sigmady2(0,jpdf2,ix2,jx2)
                              lum_gq = lum_gq
     1                               + factorQG(0,ipdf1)
     2                               * sigmady1(0,jpdf1,ix1,jx1)
     3                               * sigmady2(ipdf1,jpdf2,ix2,jx2)
                           enddo
                        endif
*     Combination and normalisation
                        sigmafkdy(jx1,jx2,jpdf1,jpdf2) =
     1                  sigmafkdy(jx1,jx2,jpdf1,jpdf2) + norm * ( 
     2                       ( cDY_NS(0,ix2,ix1,idat) 
     3                       + as * cDY_NS(1,ix2,ix1,idat) ) * lum_ns +
     4                       as * cDY_QG(1,ix2,ix1,idat) * lum_qg +
     5                       as * cDY_GQ(1,ix2,ix1,idat) * lum_gq )
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
*
*     Compute predictions
*
      if(ComputePerdictions)then
*
*     Start from the predictions obtained directly from the
*     matric elements.
*
*     Define PDFs
*
         do ix=1,nxDY
*     Projectile PDFs
            do ipdf=-6,6
               f1Q(ipdf,ix) = xPDF(ipdf,xxDY(ix)) / xxDY(ix)
            enddo
*     Target PDFs
            f2Q(0,ix)  = f1Q(0,ix)
            f2Q(1,ix)  = zarat * f1Q(ich,ix)
     1                 + ( 1d0 - zarat ) * f1Q(ich*2,ix)
            f2Q(2,ix)  = ( 1d0 - zarat ) * f1Q(ich,ix)
     1                 + zarat * f1Q(ich*2,ix)
            f2Q(-1,ix) = zarat * f1Q(-ich,ix)
     1                 + ( 1d0 - zarat ) * f1Q(-ich*2,ix)
            f2Q(-2,ix) = ( 1d0 - zarat ) * f1Q(-ich,ix)
     1                 + zarat * f1Q(-ich*2,ix)
            do ipdf=3,6
               f2Q(ipdf,ix)  = f1Q(ich*ipdf,ix)
               f2Q(-ipdf,ix) = f1Q(-ich*ipdf,ix)
            enddo
         enddo
*
*     Convolute PDFs with hard cross sections
*
         pred_me = 0d0
*     Loop on the Hard cross section grid points
         do ix1=ixp(1),nxDY
            do ix2=ixp(2),nxDY
*     Initialize luminosities
               lum_ns = 0d0
               lum_qg = 0d0
               lum_gq = 0d0
*     W^+/W^- observables
               if(obslbl(1:9).eq."EWK_WASYM".or.
     1            obslbl(1:5).eq."DYP_C")then
                  do ipdf1=-nf,nf
                     do ipdf2=-nf,nf
*     quark-antiquark luminosity
                        lum_ns = lum_ns
     1                       + factorNS(ipdf1,ipdf2)
     2                       * f1Q(ipdf1,ix1)
     3                       * f2Q(-ipdf2,ix2)
*     quark-gluon luminosity
                        lum_qg = lum_qg
     1                       + factorQG(ipdf1,ipdf2)
     2                       * f1Q(ipdf1,ix1)
     3                       * f2Q(0,ix2)    
*     gluon-quark luminosity
                        lum_gq = lum_gq
     1                       + factorQG(ipdf1,ipdf2)
     2                       * f1Q(0,ix1)
     3                       * f2Q(ipdf2,ix2)
                     enddo
                  enddo
*     gamma/Z observables
               else
                  do ipdf1=-nf,nf
                     lum_ns = lum_ns
     1                    + factorNS(ipdf1,ipdf1)
     2                    * f1Q(ipdf1,ix1)
     3                    * f2Q(-ipdf1,ix2)
                     lum_qg = lum_qg
     1                    + factorQG(ipdf1,0)
     2                    * f1Q(ipdf1,ix1)
     3                    * f2Q(0,ix2)
                     lum_gq = lum_gq
     1                    + factorQG(0,ipdf1)
     2                    * f1Q(0,ix1)
     3                    * f2Q(ipdf1,ix2)
                  enddo
               endif
*     Combination and normalisation
               pred_me = pred_me + norm * ( 
     1              ( cDY_NS(0,ix2,ix1,idat) +
     2              as * cDY_NS(1,ix2,ix1,idat) ) * lum_ns +
     3              as * cDY_QG(1,ix2,ix1,idat) * lum_qg +
     4              as * cDY_GQ(1,ix2,ix1,idat) * lum_gq )
            enddo
         enddo
*
*     Now compute the predictions using the FK tables.
*
*     Get PDFs at the initial scale
*
         call EvolveAPFEL(Q0,Q0)
         do jx=0,nx
            call xPDFall(xg(jx),xfph)
            call PDFphys2ev(xfph,xfev)
            do jpdf=1,13
               xfQ0(jpdf,jx) = xfev(jpdf)
            enddo
         enddo
*
*     Convolute with the FK tables
*
*     find ixp(1) and ixp(2)
*
         x0(1) = x1dat(idat)
         x0(2) = x2dat(idat)
         do i=1,2
            ixpfk(i) = 0
            dist(0)  = x0(i) - xg(0)
            do kx=1,nx
               dist(kx) = x0(i) - xg(kx)
               if(dist(kx)*dist(kx-1).lt.0d0) goto 101
            enddo
 101        ixpfk(i) = kx - 1
         enddo
         pred_fk = 0
         do jpdf1=1,13
            do jpdf2=1,13
               do jx1=ixpfk(1)-1,nx-1
                  do jx2=ixpfk(2)-1,nx-1
                     pred_fk = pred_fk
     1                    + sigmafkdy(jx1,jx2,jpdf1,jpdf2) 
     2                    * xfQ0(jpdf1,jx1) * xfQ0(jpdf2,jx2)
                  enddo
               enddo
            enddo
         enddo
         rel_diff = 100d0 * ( pred_me - pred_fk ) / pred_me
         write(6,"(3(a,es24.12),a)") " ME prediction: ",pred_me,
     1                               " FK prediction: ",pred_fk,
     2                               " Relative difference: ",rel_diff,
     3                               " %"
      endif
*
*     Check that the FK prediction is much bigger than the cutoff
*     used to determine the flavour map.
*
      if(dabs(cutoff/pred_fk).gt.1d-10)then
         write(6,*) "In src/FTDY/src/sigmafk_dy.f:"
         write(6,*) "the prediction is not much larger than the cutoff."
         write(6,*) "FK prediction =",pred_fk
         write(6,*) "cutoff =",cutoff
         write(6,*) "Reduce the cutoff in src/FTDY/commons/cutoff.h."
         call exit(-10)
      endif
*
      return
      end
