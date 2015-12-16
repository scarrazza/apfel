************************************************************************
*
*     EvolveQCD.f:
*
*     This routine evolves the input PDF in QCD with nf active flavours 
*     in the QCD evolution basis and replaces the input array with the
*     output one.
*
************************************************************************
      subroutine EvolveQCD(nf,fevQCD)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/EvolutionMatrices.h"
      include "../commons/MaxFlavourPDFs.h"
**
*     Input Variables
*
      integer nf
**
*     Internal Variables
*
      integer jnf
      integer i,j
      integer alpha,beta
      double precision fevQCDb(0:13,0:nint_max)
**
*     Input and Output Variables
*
      double precision fevQCD(0:13,0:nint_max)
*
*     Limit the number of flavour PDFs
*
      jnf = min(nf,nfMaxPDFs)
*
*     Apply matching conditions for the backward evolution
*
      if(sgn.eq.-1.and.nf.lt.nfi.and.nf.lt.nfMaxPDFs)then
         call MatchPDFs(nf+1,sgn,fevQCD)
      endif
*
      do alpha=0,nin(igrid)
*     Initialize backup PDF
         do i=1,13
            fevQCDb(i,alpha) = 0d0
         enddo
         do beta=0,nin(igrid)
*      Singlet and Gluon
            do i=1,2
               do j=1,2
                  fevQCDb(i,alpha) = fevQCDb(i,alpha) 
     1            + MQCDsg(nf,i,j,alpha,beta) * fevQCD(j,beta)
               enddo
            enddo
*     Total Valence
            fevQCDb(3,alpha) = fevQCDb(3,alpha) 
     1      + MQCDnsv(nf,alpha,beta) * fevQCD(3,beta)
*     V3
            fevQCDb(4,alpha) = fevQCDb(4,alpha) 
     1      + MQCDnsm(nf,alpha,beta) * fevQCD(4,beta)
*     V8
            fevQCDb(5,alpha) = fevQCDb(5,alpha) 
     1      + MQCDnsm(nf,alpha,beta) * fevQCD(5,beta)
*     V15
            if(jnf.lt.4)then
               fevQCDb(6,alpha) = fevQCDb(6,alpha) 
     1         + MQCDnsv(nf,alpha,beta) * fevQCD(3,beta)
            else
               fevQCDb(6,alpha) = fevQCDb(6,alpha) 
     1         + MQCDnsm(nf,alpha,beta) * fevQCD(6,beta)
            endif
*     V24
            if(jnf.lt.5)then
               fevQCDb(7,alpha) = fevQCDb(7,alpha) 
     1         + MQCDnsv(nf,alpha,beta) * fevQCD(3,beta)
            else
               fevQCDb(7,alpha) = fevQCDb(7,alpha) 
     1         + MQCDnsm(nf,alpha,beta) * fevQCD(7,beta)
            endif
*     V35
            if(jnf.lt.6)then
               fevQCDb(8,alpha) = fevQCDb(8,alpha) 
     1         + MQCDnsv(nf,alpha,beta) * fevQCD(3,beta)
            else
               fevQCDb(8,alpha) = fevQCDb(8,alpha) 
     1         + MQCDnsm(nf,alpha,beta) * fevQCD(8,beta)
            endif
*     T3
            fevQCDb(9,alpha) = fevQCDb(9,alpha) 
     1      + MQCDnsp(nf,alpha,beta) * fevQCD(9,beta)
*     T8
            fevQCDb(10,alpha) = fevQCDb(10,alpha) 
     1      + MQCDnsp(nf,alpha,beta) * fevQCD(10,beta)
*     T15
            if(jnf.lt.4)then
               do j=1,2
                  fevQCDb(11,alpha) = fevQCDb(11,alpha) 
     1            + MQCDsg(nf,1,j,alpha,beta) * fevQCD(j,beta)
               enddo
            else
               fevQCDb(11,alpha) = fevQCDb(11,alpha) 
     1         + MQCDnsp(nf,alpha,beta) * fevQCD(11,beta)
            endif
*     T24
            if(jnf.lt.5)then
               do j=1,2
                  fevQCDb(12,alpha) = fevQCDb(12,alpha) 
     1            + MQCDsg(nf,1,j,alpha,beta) * fevQCD(j,beta)
               enddo
            else
               fevQCDb(12,alpha) = fevQCDb(12,alpha) 
     1         + MQCDnsp(nf,alpha,beta) * fevQCD(12,beta)
            endif
*     T35
            if(jnf.lt.6)then
               do j=1,2
                  fevQCDb(13,alpha) = fevQCDb(13,alpha) 
     1            + MQCDsg(nf,1,j,alpha,beta) * fevQCD(j,beta)
               enddo
            else
               fevQCDb(13,alpha) = fevQCDb(13,alpha) 
     1         + MQCDnsp(nf,alpha,beta) * fevQCD(13,beta)
            endif
         enddo
      enddo
*
*     Apply matching conditions for the forward evolution
*
      if(sgn.eq.1.and.nf.lt.nff.and.nf.lt.nfMaxPDFs)then
         call MatchPDFs(nf+1,sgn,fevQCDb)
      endif
*
*     Copy backup PDFs into main PDFs
*
      do alpha=0,nin(igrid)
         do i=1,13
            fevQCD(i,alpha) = fevQCDb(i,alpha)
         enddo
      enddo
*
      return
      end
