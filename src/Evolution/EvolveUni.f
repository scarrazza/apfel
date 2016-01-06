************************************************************************
*
*     EvolveUni.f:
*
*     This routine evolves the input PDF in QCD and QED with nf active
*     flavours in the unified evolution basis and replaces the input array
*     with the output one.
*
************************************************************************
      subroutine EvolveUni(nf,nl,flevUni,fevUni)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/EvolutionMatrices.h"
      include "../commons/MaxFlavourPDFs.h"
**
*     Input Variables
*
      integer nf,nl
**
*     Internal Variables
*
      integer jnf
      integer i,j
      integer alpha,beta
      double precision fevUnib(0:13,0:nint_max)
      double precision flevUnib(6,0:nint_max)
      double precision fevQCD(0:13,0:nint_max)
      double precision fsg1(5)
**
*     Input and Output Variables
*
      double precision flevUni(6,0:nint_max)
      double precision fevUni(0:13,0:nint_max)
*
*     Limit the number of flavour PDFs
*
      jnf = min(nf,nfMaxPDFs)
*
*     Apply matching conditions for the backward evolution
*
      if(sgn.eq.-1.and.nf.lt.nfli(nl).and.nf.lt.nfMaxPDFs)then
*     Rotate to the QCD evolution basis
         call PDFevUni2evQCD(fevUni,fevQCD)
*     Apply matching conditions
         call MatchPDFs(nf+1,sgn,fevQCD)
*     Rotate back to the unified evolution basis
         call PDFevQCD2evUni(fevQCD,fevUni)
      endif
*
      do alpha=0,nin(igrid)
*     Initialize backup PDFs
         do i=0,13
            fevUnib(i,alpha) = 0d0
         enddo
         do i=1,6
            flevUnib(i,alpha) = 0d0
         enddo
         do beta=0,nin(igrid)
*      Singlet 1
            fsg1(1) = fevUni(0,beta)
            fsg1(2) = fevUni(1,beta)
            fsg1(3) = fevUni(2,beta)
            fsg1(4) = fevUni(3,beta)
            fsg1(5) = flevUni(1,beta)
            do i=1,4
               do j=1,5
                  fevUnib(i-1,alpha) = fevUnib(i-1,alpha) 
     1            + MUnisg1(nf,nl,i,j,alpha,beta) * fsg1(j)
               enddo
            enddo
            do j=1,5
               flevUnib(1,alpha) = flevUnib(1,alpha) 
     1         + MUnisg1(nf,nl,5,j,alpha,beta) * fsg1(j)
            enddo
*      Singlet 2
            do i=1,2
               do j=1,2
                  fevUnib(i+7,alpha) = fevUnib(i+7,alpha) 
     1            + MUnisg2(nf,nl,i,j,alpha,beta) * fevUni(j+7,beta)
               enddo
            enddo
*     Tu1
            if(jnf.ge.4)then
               fevUnib(4,alpha) = fevUnib(4,alpha)
     1         + MUninspu(nf,nl,alpha,beta) * fevUni(4,beta)
            else
               do i=1,5
                  fevUnib(4,alpha) = fevUnib(4,alpha)
     1                             + ( MUnisg1(nf,nl,3,i,alpha,beta)
     2                             + MUnisg1(nf,nl,4,i,alpha,beta) )
     3                             * fsg1(i) / 2d0
               enddo
            endif
*     Tu2
            if(jnf.ge.6)then
               fevUnib(5,alpha) = fevUnib(5,alpha)
     1         + MUninspu(nf,nl,alpha,beta) * fevUni(5,beta)
            else
               do i=1,5
                  fevUnib(5,alpha) = fevUnib(5,alpha)
     1                             + ( MUnisg1(nf,nl,3,i,alpha,beta)
     2                             + MUnisg1(nf,nl,4,i,alpha,beta) )
     3                             * fsg1(i) / 2d0
               enddo
            endif
*     Td1
            fevUnib(6,alpha) = fevUnib(6,alpha)
     1      + MUninspd(nf,nl,alpha,beta) * fevUni(6,beta)
*     Td2
            if(jnf.ge.5)then
               fevUnib(7,alpha) = fevUnib(7,alpha)
     1         + MUninspd(nf,nl,alpha,beta) * fevUni(7,beta)
            else
               do i=1,5
                  fevUnib(7,alpha) = fevUnib(7,alpha)
     1                             + ( MUnisg1(nf,nl,3,i,alpha,beta)
     2                             - MUnisg1(nf,nl,4,i,alpha,beta) )
     3                             * fsg1(i) / 2d0
               enddo
            endif
*     Vu1
            if(jnf.ge.4)then
               fevUnib(10,alpha) = fevUnib(10,alpha)
     1         + MUninsmu(nf,nl,alpha,beta) * fevUni(10,beta)
            else
               do i=1,2
                  fevUnib(10,alpha) = fevUnib(10,alpha)
     1                              + ( MUnisg2(nf,nl,1,i,alpha,beta)
     2                              + MUnisg2(nf,nl,2,i,alpha,beta) )
     3                              * fevUni(i+7,beta) / 2d0
               enddo
            endif
*     Vu2
            if(jnf.ge.6)then
               fevUnib(11,alpha) = fevUnib(11,alpha)
     1         + MUninsmu(nf,nl,alpha,beta) * fevUni(11,beta)
            else
               do i=1,2
                  fevUnib(11,alpha) = fevUnib(11,alpha)
     1                              + ( MUnisg2(nf,nl,1,i,alpha,beta)
     2                              + MUnisg2(nf,nl,2,i,alpha,beta) )
     3                              * fevUni(i+7,beta) / 2d0
               enddo
            endif
*     Vd1
            fevUnib(12,alpha) = fevUnib(12,alpha)
     1      + MUninsmd(nf,nl,alpha,beta) * fevUni(12,beta)
*     Vd2
            if(jnf.ge.5)then
               fevUnib(13,alpha) = fevUnib(13,alpha)
     1         + MUninsmd(nf,nl,alpha,beta) * fevUni(13,beta)
            else
               do i=1,2
                  fevUnib(13,alpha) = fevUnib(13,alpha)
     1                              + ( MUnisg2(nf,nl,1,i,alpha,beta)
     2                              - MUnisg2(nf,nl,2,i,alpha,beta) )
     3                              * fevUni(i+7,beta) / 2d0
               enddo
            endif
*     Tl3
            flevUnib(2,alpha) = flevUnib(2,alpha)
     1      + MUninslep(nl,alpha,beta) * flevUni(2,beta)
*     Tl8
            if(nl.ge.3)then
               flevUnib(3,alpha) = flevUnib(3,alpha)
     1         + MUninslep(nl,alpha,beta) * flevUni(3,beta)
            else
               do j=1,5
                  flevUnib(3,alpha) = flevUnib(3,alpha) 
     1            + MUnisg1(nf,nl,5,j,alpha,beta) * fsg1(j)
               enddo
            endif
*     Vl
            flevUnib(4,alpha) = flevUnib(4,alpha)
     1      + MUninslep(nl,alpha,beta) * flevUni(5,beta)
*     Vl3
            flevUnib(5,alpha) = flevUnib(5,alpha)
     1      + MUninslep(nl,alpha,beta) * flevUni(5,beta)
*     Vl8
            flevUnib(6,alpha) = flevUnib(6,alpha)
     1      + MUninslep(nl,alpha,beta) * flevUni(6,beta)
         enddo
      enddo
*
*     Apply matching conditions for the forward evolution
*
      if(sgn.eq.1.and.nf.lt.nflf(nl).and.nf.lt.nfMaxPDFs)then
*     Rotate to the QCD evolution basis
         call PDFevUni2evQCD(fevUnib,fevQCD)
*     Apply matching conditions
         call MatchPDFs(nf+1,sgn,fevQCD)
*     Rotate back to the unified evolution basis
         call PDFevQCD2evUni(fevQCD,fevUnib)
      endif
*
*     Copy backup PDFs into main PDFs
*
      do alpha=0,nin(igrid)
         do i=0,13
            fevUni(i,alpha) = fevUnib(i,alpha)
         enddo
         do i=1,6
            flevUni(i,alpha) = flevUnib(i,alpha)
         enddo
      enddo
*
      return
      end
