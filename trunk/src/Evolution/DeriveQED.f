************************************************************************
*
*     DeriveQED.f:
*
*     This routine evolves the input PDF in QED in the QED evolution
*     basis and replaces the input array with the output one.
*
************************************************************************
      subroutine DeriveQED(dfevQED)
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
      double precision dfevQEDb(0:13,0:nint_max)
**
*     Input and Output Variables
*
      double precision dfevQED(0:13,0:nint_max)
*
      do alpha=0,nin(igrid)
*     Initialize backup PDF
         do i=1,13
            dfevQEDb(i,alpha) = 0d0
         enddo
         do beta=0,nin(igrid)
*     Gamma, Singlet and Detla
            do i=1,3
               do j=1,3
                  dfevQEDb(i,alpha) = dfevQEDb(i,alpha) 
     1            +  dMQEDsg(i,j,alpha,beta) * dfevQED(j,beta)
               enddo
            enddo
*     Duc
            if(wnf.lt.4)then
               do j=1,3
                  dfevQEDb(4,alpha) = dfevQEDb(4,alpha) 
     1            +  ( dMQEDsg(2,j,alpha,beta) 
     2            +    dMQEDsg(3,j,alpha,beta) ) 
     3            * dfevQED(j,beta) / 2d0
               enddo
            else
               dfevQEDb(4,alpha) = dfevQEDb(4,alpha) 
     1         +  dMQEDnsp(alpha,beta) * dfevQED(4,beta)
            endif
*     Dds
            dfevQEDb(5,alpha) = dfevQEDb(5,alpha) 
     1      +  dMQEDnsm(alpha,beta) * dfevQED(5,beta)
*     Dsb
            if(wnf.lt.5)then
               do j=1,3
                  dfevQEDb(6,alpha) = dfevQEDb(6,alpha)
     1            +  ( dMQEDsg(2,j,alpha,beta)
     2            -    dMQEDsg(3,j,alpha,beta) )
     3            * dfevQED(j,beta) / 4d0
               enddo
               dfevQEDb(6,alpha) = dfevQEDb(6,alpha) 
     1         - dMQEDnsm(alpha,beta) * dfevQED(5,beta) / 2d0
            else
               dfevQEDb(6,alpha) = dfevQEDb(6,alpha) 
     1         +  dMQEDnsm(alpha,beta) * dfevQED(6,beta)
            endif
*     Dct
            if(wnf.lt.4)then
               dfevQEDb(7,alpha) = 0d0
            elseif(wnf.lt.6)then
               do j=1,3
                  dfevQEDb(7,alpha) = dfevQEDb(7,alpha)
     1            +  ( dMQEDsg(2,j,alpha,beta)
     2            +    dMQEDsg(3,j,alpha,beta) )
     3            * dfevQED(j,beta) / 4d0
               enddo
               dfevQEDb(7,alpha) = dfevQEDb(7,alpha) 
     1         - dMQEDnsp(alpha,beta) * dfevQED(4,beta) / 2d0
            else
               dfevQEDb(7,alpha) = dfevQEDb(7,alpha) 
     1         +  dMQEDnsp(alpha,beta) * dfevQED(7,beta)
            endif
*     u-
            dfevQEDb(8,alpha) = dfevQEDb(8,alpha) 
     1      +  dMQEDnsp(alpha,beta) * dfevQED(8,beta)
*     d-
            dfevQEDb(9,alpha) = dfevQEDb(9,alpha) 
     1      +  dMQEDnsm(alpha,beta) * dfevQED(9,beta)
*     s-
             dfevQEDb(10,alpha) = dfevQEDb(10,alpha) 
     1       +  dMQEDnsm(alpha,beta) * dfevQED(10,beta)
*     c-
            if(wnf.lt.4)then
               dfevQEDb(11,alpha) = 0d0
            else
               dfevQEDb(11,alpha) = dfevQEDb(11,alpha) 
     1         +  dMQEDnsp(alpha,beta) * dfevQED(11,beta)
            endif
*     b-
            if(wnf.lt.5)then
               dfevQEDb(12,alpha) = 0d0
            else
               dfevQEDb(12,alpha) = dfevQEDb(12,alpha) 
     1         +  dMQEDnsm(alpha,beta) * dfevQED(12,beta)
            endif
*     t-
            if(wnf.lt.6)then
               dfevQEDb(13,alpha) = 0d0
            else
               dfevQEDb(13,alpha) = dfevQEDb(13,alpha) 
     1         +  dMQEDnsp(alpha,beta) * dfevQED(13,beta)
            endif
         enddo
      enddo
*
*     Copy backup PDFs into main PDFs
*
      do alpha=0,nin(igrid)
         do i=1,13
            dfevQED(i,alpha)  = dfevQEDb(i,alpha)
         enddo
      enddo
*
      return
      end
