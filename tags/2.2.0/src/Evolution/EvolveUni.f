************************************************************************
*
*     EvolveUni.f:
*
*     This routine evolves the input PDF in QCD and QED with nf active
*     flavours in the unified evolution basis and replaces the input array
*     with the output one.
*
************************************************************************
      subroutine EvolveUni(nf,fevUni)
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
      double precision fevUnib(0:13,0:nint_max)
      double precision fevQCD(0:13,0:nint_max)
**
*     Input and Output Variables
*
      double precision fevUni(0:13,0:nint_max)
*
*     Apply matching conditions for the backward evolution
*
      if(sgn.eq.-1)then
         if(nf.lt.nfi)then
*     Rotate to the QCD evolution basis
            call PDFevUni2evQCD(fevUni,fevQCD)
*     Apply matching conditions
            call MatchPDFs(nf,fevQCD)
*     Rotate back to the unified evolution basis
            call PDFevQCD2evUni(fevQCD,fevUni)
         endif
      endif
*
      do alpha=0,nin(igrid)
*     Initialize backup PDFs
         do i=0,13
            fevUnib(i,alpha) = 0d0
         enddo
         do beta=0,nin(igrid)
*      Singlet 1
            do i=1,4
               do j=1,4
                  fevUnib(i-1,alpha) = fevUnib(i-1,alpha) 
     1            + MUnisg1(nf,i,j,alpha,beta) * fevUni(j-1,beta)
               enddo
            enddo
*      Singlet 2
            do i=1,2
               do j=1,2
                  fevUnib(i+7,alpha) = fevUnib(i+7,alpha) 
     1            + MUnisg2(nf,i,j,alpha,beta) * fevUni(j+7,beta)
               enddo
            enddo
*     Tu1
            if(nf.ge.4)then
               fevUnib(4,alpha) = fevUnib(4,alpha)
     1         + MUninspu(nf,alpha,beta) * fevUni(4,beta)
            else
               do i=1,4
                  fevUnib(4,alpha) = fevUnib(4,alpha)
     1                             + ( MUnisg1(nf,3,i,alpha,beta)
     2                             + MUnisg1(nf,4,i,alpha,beta) )
     3                             * fevUni(i-1,beta) / 2d0
               enddo
            endif
*     Tu2
            if(nf.ge.6)then
               fevUnib(5,alpha) = fevUnib(5,alpha)
     1         + MUninspu(nf,alpha,beta) * fevUni(5,beta)
            else
               do i=1,4
                  fevUnib(5,alpha) = fevUnib(5,alpha)
     1                             + ( MUnisg1(nf,3,i,alpha,beta)
     2                             + MUnisg1(nf,4,i,alpha,beta) )
     3                             * fevUni(i-1,beta) / 2d0
               enddo
            endif
*     Td1
            fevUnib(6,alpha) = fevUnib(6,alpha)
     1      + MUninspd(nf,alpha,beta) * fevUni(6,beta)
*     Td2
            if(nf.ge.5)then
               fevUnib(7,alpha) = fevUnib(7,alpha)
     1         + MUninspd(nf,alpha,beta) * fevUni(7,beta)
            else
               do i=1,4
                  fevUnib(7,alpha) = fevUnib(7,alpha)
     1                             + ( MUnisg1(nf,3,i,alpha,beta)
     2                             - MUnisg1(nf,4,i,alpha,beta) )
     3                             * fevUni(i-1,beta) / 2d0
               enddo
            endif
*     Vu1
            if(nf.ge.4)then
               fevUnib(10,alpha) = fevUnib(10,alpha)
     1         + MUninsmu(nf,alpha,beta) * fevUni(10,beta)
            else
               do i=1,2
                  fevUnib(10,alpha) = fevUnib(10,alpha)
     1                              + ( MUnisg2(nf,1,i,alpha,beta)
     2                              + MUnisg2(nf,2,i,alpha,beta) )
     3                              * fevUni(i+7,beta) / 2d0
               enddo
            endif
*     Vu2
            if(nf.ge.6)then
               fevUnib(11,alpha) = fevUnib(11,alpha)
     1         + MUninsmu(nf,alpha,beta) * fevUni(11,beta)
            else
               do i=1,2
                  fevUnib(11,alpha) = fevUnib(11,alpha)
     1                              + ( MUnisg2(nf,1,i,alpha,beta)
     2                              + MUnisg2(nf,2,i,alpha,beta) )
     3                              * fevUni(i+7,beta) / 2d0
               enddo
            endif
*     Vd1
            fevUnib(12,alpha) = fevUnib(12,alpha)
     1      + MUninsmd(nf,alpha,beta) * fevUni(12,beta)
*     Vd2
            if(nf.ge.5)then
               fevUnib(13,alpha) = fevUnib(13,alpha)
     1         + MUninsmd(nf,alpha,beta) * fevUni(13,beta)
            else
               do i=1,2
                  fevUnib(13,alpha) = fevUnib(13,alpha)
     1                              + ( MUnisg2(nf,1,i,alpha,beta)
     2                              - MUnisg2(nf,2,i,alpha,beta) )
     3                              * fevUni(i+7,beta) / 2d0
               enddo
            endif
         enddo
      enddo
*
*     Apply matching conditions for the forward evolution
*
      if(sgn.eq.1)then
         if(nf.lt.nff)then
*     Rotate to the QCD evolution basis
            call PDFevUni2evQCD(fevUnib,fevQCD)
*     Apply matching conditions
            call MatchPDFs(nf+1,fevQCD)
*     Rotate back to the unified evolution basis
            call PDFevQCD2evUni(fevQCD,fevUnib)
         endif
      endif
*
*     Copy backup PDFs into main PDFs
*
      do alpha=0,nin(igrid)
         do i=0,13
            fevUni(i,alpha) = fevUnib(i,alpha)
         enddo
      enddo
*
      return
      end
