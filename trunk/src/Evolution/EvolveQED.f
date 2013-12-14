************************************************************************
*
*     EvolveQED.f:
*
*     This routine evolves the input PDF in QED with nf active flavours 
*     in the QED evolution basis and replaces the input array with the
*     output one.
*
************************************************************************
      subroutine EvolveQED(nf,fevQED)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/EvolutionMatrices.h"
**
*     Input Variables
*
      integer nf
**
*     Internal Variables
*
      integer i,j
      integer alpha,beta
      double precision fevQEDb(0:13,0:nint_max)
**
*     Input and Output Variables
*
      double precision fevQED(0:13,0:nint_max)
*
      do alpha=0,nin(igrid)
*     Initialize backup PDF
         do i=1,13
            fevQEDb(i,alpha) = 0d0
         enddo
         do beta=0,nin(igrid)
*     Gamma, Singlet and Detla
            do i=1,3
               do j=1,3
                  fevQEDb(i,alpha) = fevQEDb(i,alpha) 
     1            +  MQEDsg(nf,i,j,alpha,beta) * fevQED(j,beta)
               enddo
            enddo
*     Duc
            if(nf.lt.4)then
               do j=1,3
                  fevQEDb(4,alpha) = fevQEDb(4,alpha) 
     1            +  ( MQEDsg(nf,2,j,alpha,beta) 
     2            +    MQEDsg(nf,3,j,alpha,beta) ) 
     3            * fevQED(j,beta) / 2d0
               enddo
            else
               fevQEDb(4,alpha) = fevQEDb(4,alpha) 
     1         +  MQEDnsp(nf,alpha,beta) * fevQED(4,beta)
            endif
*     Dds
            fevQEDb(5,alpha) = fevQEDb(5,alpha) 
     1      +  MQEDnsm(nf,alpha,beta) * fevQED(5,beta)
*     Dsb
            if(nf.lt.5)then
               do j=1,3
                  fevQEDb(6,alpha) = fevQEDb(6,alpha)
     1            +  ( MQEDsg(nf,2,j,alpha,beta)
     2            -    MQEDsg(nf,3,j,alpha,beta) )
     3            * fevQED(j,beta) / 4d0
               enddo
               fevQEDb(6,alpha) = fevQEDb(6,alpha) 
     1         - MQEDnsm(nf,alpha,beta) * fevQED(5,beta) / 2d0
            else
               fevQEDb(6,alpha) = fevQEDb(6,alpha) 
     1         +  MQEDnsm(nf,alpha,beta) * fevQED(6,beta)
            endif
*     Dct
            if(nf.lt.4)then
               fevQEDb(7,alpha) = 0d0
            elseif(nf.lt.6)then
               do j=1,3
                  fevQEDb(7,alpha) = fevQEDb(7,alpha)
     1            +  ( MQEDsg(nf,2,j,alpha,beta)
     2            +    MQEDsg(nf,3,j,alpha,beta) )
     3            * fevQED(j,beta) / 4d0
               enddo
               fevQEDb(7,alpha) = fevQEDb(7,alpha) 
     1         - MQEDnsp(nf,alpha,beta) * fevQED(4,beta) / 2d0
            else
               fevQEDb(7,alpha) = fevQEDb(7,alpha) 
     1         +  MQEDnsp(nf,alpha,beta) * fevQED(7,beta)
            endif
*     u-
            fevQEDb(8,alpha) = fevQEDb(8,alpha) 
     1      +  MQEDnsp(nf,alpha,beta) * fevQED(8,beta)
*     d-
            fevQEDb(9,alpha) = fevQEDb(9,alpha) 
     1      +  MQEDnsm(nf,alpha,beta) * fevQED(9,beta)
*     s-
             fevQEDb(10,alpha) = fevQEDb(10,alpha) 
     1       +  MQEDnsm(nf,alpha,beta) * fevQED(10,beta)
*     c-
            if(nf.lt.4)then
               fevQEDb(11,alpha) = 0d0
            else
               fevQEDb(11,alpha) = fevQEDb(11,alpha) 
     1         +  MQEDnsp(nf,alpha,beta) * fevQED(11,beta)
            endif
*     b-
            if(nf.lt.5)then
               fevQEDb(12,alpha) = 0d0
            else
               fevQEDb(12,alpha) = fevQEDb(12,alpha) 
     1         +  MQEDnsm(nf,alpha,beta) * fevQED(12,beta)
            endif
*     t-
            if(nf.lt.6)then
               fevQEDb(13,alpha) = 0d0
            else
               fevQEDb(13,alpha) = fevQEDb(13,alpha) 
     1         +  MQEDnsp(nf,alpha,beta) * fevQED(13,beta)
            endif
         enddo
      enddo
*
*     Copy backup PDFs into main PDFs
*
      do alpha=0,nin(igrid)
         do i=1,13
            fevQED(i,alpha)  = fevQEDb(i,alpha)
         enddo
      enddo
*
      return
      end
