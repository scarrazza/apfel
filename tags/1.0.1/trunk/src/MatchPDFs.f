************************************************************************
*
*     MatchPDFs.f:
*
*     This routine mathches PDF in QCD at the nf-th heavy quark
*     threshold in the evolution basis and replaces the input array with 
*     the output one.
*
************************************************************************
      subroutine MatchPDFs(nf,fevQCD)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf
**
*     Internal Variables
*
      integer i,j
      integer alpha,beta
      integer mapp(2,2)
      double precision coup,a_QCD
      double precision integralsMatching
      double precision integ(5,0:nint_max)
      double precision fevQCDb(0:13,0:nint_max)
**
*     Input and Output Variables
*
      double precision fevQCD(0:13,0:nint_max)
*
*     Get alphas value at the heavy quark threshold (with nf active flavours)
*
      coup = a_QCD(m2th(nf))
*
*     Contruct matching conditions at this threshod
*
      do alpha=0,nin(igrid)
         do i=1,5
            integ(i,alpha) = integralsMatching(0,alpha,coup,i)
         enddo
      enddo
*
*     Singlet map
*
      mapp(1,1) = 2
      mapp(1,2) = 3
      mapp(2,1) = 4
      mapp(2,2) = 5
*
      do alpha=0,nin(igrid)
*     Initialize backup PDF
         do i=1,13
            fevQCDb(i,alpha) = 0d0
         enddo
         do beta=0,nin(igrid)-alpha
*      Singlet and Gluon
            do i=1,2
               do j=1,2
                  fevQCDb(i,alpha) = fevQCDb(i,alpha) 
     1            + integ(mapp(i,j),beta) * fevQCD(j,alpha+beta)
               enddo
            enddo
*     Total Valence
            fevQCDb(3,alpha) = fevQCDb(3,alpha) 
     1                       + integ(1,beta) * fevQCD(3,alpha+beta)
*     V3
            fevQCDb(4,alpha) = fevQCDb(4,alpha) 
     1                       + integ(1,beta) * fevQCD(4,alpha+beta)
*     V8
            fevQCDb(5,alpha) = fevQCDb(5,alpha) 
     1                       + integ(1,beta) * fevQCD(5,alpha+beta)
*     V15
            fevQCDb(6,alpha) = fevQCDb(6,alpha) 
     1                       + integ(1,beta) * fevQCD(6,alpha+beta)
*     V24
            fevQCDb(7,alpha) = fevQCDb(7,alpha) 
     1                       + integ(1,beta) * fevQCD(7,alpha+beta)
*     V35
            fevQCDb(8,alpha) = fevQCDb(8,alpha) 
     1                       + integ(1,beta) * fevQCD(8,alpha+beta)
*     T3
            fevQCDb(9,alpha) = fevQCDb(9,alpha) 
     1                       + integ(1,beta) * fevQCD(9,alpha+beta)
*     T8
            fevQCDb(10,alpha) = fevQCDb(10,alpha) 
     1                        + integ(1,beta) * fevQCD(10,alpha+beta)
*     T15
            if(nf.eq.4)then
               fevQCDb(11,alpha) = fevQCDb(11,alpha) 
     1         + ( integ(1,beta) - 3d0 * ( integ(2,beta) 
     2         - integ(1,beta) ) )  * fevQCD(1,alpha+beta)
     3         - 3d0 * integ(3,beta) * fevQCD(2,alpha+beta)
            else
               fevQCDb(11,alpha) = fevQCDb(11,alpha) 
     1                           + integ(1,beta) * fevQCD(11,alpha+beta)
            endif
*     T24
            if(nf.eq.4)then
               fevQCDb(12,alpha) = fevQCDb(12,alpha) 
     1                           + integ(2,beta) * fevQCD(1,alpha+beta)
     2                           + integ(3,beta) * fevQCD(2,alpha+beta)
            elseif(nf.eq.5)then
               fevQCDb(12,alpha) = fevQCDb(12,alpha) 
     1         + ( integ(1,beta) - 4d0 * ( integ(2,beta) 
     2         - integ(1,beta) ) )  * fevQCD(1,alpha+beta)
     3         - 4d0 * integ(3,beta) * fevQCD(2,alpha+beta)
            else
               fevQCDb(12,alpha) = fevQCDb(12,alpha) 
     1                           + integ(1,beta) * fevQCD(12,alpha+beta)
            endif
*     T35
            if(nf.eq.6)then
               fevQCDb(13,alpha) = fevQCDb(13,alpha) 
     1         + ( integ(1,beta) - 5d0 * ( integ(2,beta) 
     2         - integ(1,beta) ) )  * fevQCD(1,alpha+beta)
     3         - 5d0 * integ(3,beta) * fevQCD(2,alpha+beta)
            else
               fevQCDb(13,alpha) = fevQCDb(13,alpha) 
     1                           + integ(2,beta) * fevQCD(1,alpha+beta)
     2                           + integ(3,beta) * fevQCD(2,alpha+beta)
            endif
         enddo
      enddo
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
