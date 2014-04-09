************************************************************************
*
*     DeriveQCD.f:
*
*     This routine derives the input PDF in QCD in the QCD evolution
*     basis and replaces the input array with the output one.
*
************************************************************************
      subroutine DeriveQCD(dfevQCD)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/EvolutionMatrices.h"
      include "../commons/wrap.h"
**
*     Internal Variables
*
      integer i,j
      integer alpha,beta
      double precision dfevQCDb(0:13,0:nint_max)
**
*     Input and Output Variables
*
      double precision dfevQCD(0:13,0:nint_max)
*
      do alpha=0,nin(igrid)
*     Initialize backup PDF
         do i=1,13
            dfevQCDb(i,alpha) = 0d0
         enddo
         do beta=0,nin(igrid)
*      Singlet and Gluon
            do i=1,2
               do j=1,2
                  dfevQCDb(i,alpha) = dfevQCDb(i,alpha) 
     1            + dMQCDsg(i,j,alpha,beta) * dfevQCD(j,beta)
               enddo
            enddo
*     Total Valence
            dfevQCDb(3,alpha) = dfevQCDb(3,alpha) 
     1      + dMQCDnsv(alpha,beta) * dfevQCD(3,beta)
*     V3
            dfevQCDb(4,alpha) = dfevQCDb(4,alpha) 
     1      + dMQCDnsm(alpha,beta) * dfevQCD(4,beta)
*     V8
            dfevQCDb(5,alpha) = dfevQCDb(5,alpha) 
     1      + dMQCDnsm(alpha,beta) * dfevQCD(5,beta)
*     V15
            if(wnf.lt.4)then
               dfevQCDb(6,alpha) = dfevQCDb(6,alpha) 
     1         + dMQCDnsv(alpha,beta) * dfevQCD(3,beta)
            else
               dfevQCDb(6,alpha) = dfevQCDb(6,alpha) 
     1         + dMQCDnsm(alpha,beta) * dfevQCD(6,beta)
            endif
*     V24
            if(wnf.lt.5)then
               dfevQCDb(7,alpha) = dfevQCDb(7,alpha) 
     1         + dMQCDnsv(alpha,beta) * dfevQCD(3,beta)
            else
               dfevQCDb(7,alpha) = dfevQCDb(7,alpha) 
     1         + dMQCDnsm(alpha,beta) * dfevQCD(7,beta)
            endif
*     V35
            if(wnf.lt.6)then
               dfevQCDb(8,alpha) = dfevQCDb(8,alpha) 
     1         + dMQCDnsv(alpha,beta) * dfevQCD(3,beta)
            else
               dfevQCDb(8,alpha) = dfevQCDb(8,alpha) 
     1         + dMQCDnsm(alpha,beta) * dfevQCD(8,beta)
            endif
*     T3
            dfevQCDb(9,alpha) = dfevQCDb(9,alpha) 
     1      + dMQCDnsp(alpha,beta) * dfevQCD(9,beta)
*     T8
            dfevQCDb(10,alpha) = dfevQCDb(10,alpha) 
     1      + dMQCDnsp(alpha,beta) * dfevQCD(10,beta)
*     T15
            if(wnf.lt.4)then
               do j=1,2
                  dfevQCDb(11,alpha) = dfevQCDb(11,alpha) 
     1            + dMQCDsg(1,j,alpha,beta) * dfevQCD(j,beta)
               enddo
            else
               dfevQCDb(11,alpha) = dfevQCDb(11,alpha) 
     1         + dMQCDnsp(alpha,beta) * dfevQCD(11,beta)
            endif
*     T24
            if(wnf.lt.5)then
               do j=1,2
                  dfevQCDb(12,alpha) = dfevQCDb(12,alpha) 
     1            + dMQCDsg(1,j,alpha,beta) * dfevQCD(j,beta)
               enddo
            else
               dfevQCDb(12,alpha) = dfevQCDb(12,alpha) 
     1         + dMQCDnsp(alpha,beta) * dfevQCD(12,beta)
            endif
*     T35
            if(wnf.lt.6)then
               do j=1,2
                  dfevQCDb(13,alpha) = dfevQCDb(13,alpha) 
     1            + dMQCDsg(1,j,alpha,beta) * dfevQCD(j,beta)
               enddo
            else
               dfevQCDb(13,alpha) = dfevQCDb(13,alpha) 
     1         + dMQCDnsp(alpha,beta) * dfevQCD(13,beta)
            endif
         enddo
      enddo
*
*     Copy backup PDFs into main PDFs
*
      do alpha=0,nin(igrid)
         do i=1,13
            dfevQCD(i,alpha) = dfevQCDb(i,alpha)
         enddo
      enddo
*
      return
      end
