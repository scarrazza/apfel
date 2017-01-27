************************************************************************
*
*     ExternalEvolutionOperator.f:
*
*     This function returns the "external" evolution operator on the
*     interpolation grid. This is possible only if the internal grids are 
*     used or if one single external grid is provided.
*
************************************************************************
      function ExternalEvolutionOperator(Bs2Bs,i,j,x,beta)
*
      implicit none
*
      include "../commons/EvolOp.h"
      include "../commons/grid.h"
      include "../commons/EvolutionOperator.h"
**
*     Input Variables
*
      integer i,j
      integer beta
      double precision x
      character*5 Bs2Bs
**
*     Internal Variables
*
      integer n
      integer alpha
      integer lx,ux
      double precision w_int_gen,wg
      double precision tol
      parameter(tol=1d-10)
**
*     Output Variables
*
      double precision ExternalEvolutionOperator
*
*     Check whether the evolution operator has been actually computed
*
      if(.not.EvolOp)then
         write(6,*) "The evolution operator computation is disabled."
         write(6,*) "The 'ExternalEvolutionOperator' function cannot",
     1              " be used."
         write(6,*) "   "
         call exit(-10)
      endif
*
*     Check consistency of the input variables
*
      if(Bs2Bs.ne."Ev2Ev".and.
     1   Bs2Bs.ne."Ev2Ph".and.
     2   Bs2Bs.ne."Ph2Ph")then
         write(6,*) "In ExternalEvolutionOperator.f:"
         write(6,*) "Invalid Basis flag, Bs2Bs = ",Bs2Bs
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'Ev2Ev'"
         write(6,*) "- 'Ev2Ph'"
         write(6,*) "- 'Ph2Ph'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      if(Bs2Bs.eq."Ev2Ev")then
         if(i.lt.0.or.i.gt.13)then
            write(6,*) "In ExternalEvolutionOperator.f:"
            write(6,*) "Invalid index, i =",i
            call exit(-10)
         endif
         if(j.lt.0.or.j.gt.13)then
            write(6,*) "In ExternalEvolutionOperator.f:"
            write(6,*) "Invalid index, j =",j
            call exit(-10)
         endif
      elseif(Bs2Bs.eq."Ev2Ph")then
         if(i.lt.-7.or.i.gt.6)then
            write(6,*) "In ExternalEvolutionOperator.f:"
            write(6,*) "Invalid index, i =",i
            call exit(-10)
         endif
         if(j.lt.0.or.j.gt.13)then
            write(6,*) "In ExternalEvolutionOperator.f:"
            write(6,*) "Invalid index, j =",j
            call exit(-10)
         endif
      elseif(Bs2Bs.eq."Ph2Ph")then
         if(i.lt.-7.or.i.gt.6)then
            write(6,*) "In ExternalEvolutionOperator.f:"
            write(6,*) "Invalid index, i =",i
            call exit(-10)
         endif
         if(j.lt.-7.or.j.gt.6)then
            write(6,*) "In ExternalEvolutionOperator.f:"
            write(6,*) "Invalid index, j =",j
            call exit(-10)
         endif
      endif
*
      if(x.lt.xmin(1)-tol.or.x.gt.xmax+tol)then
         write(6,*) "In ExternalEvolutionOperator.f:"
         write(6,*) "Invalid value of x =",x
         call exit(-10)
      endif
      if (x.lt.xmin(1)) x = xmin(1)
      if (x.gt.xmax) x = 1d0
      if(beta.lt.0.or.beta.gt.nin(0))then
         write(6,*) "In ExternalEvolutionOperator.f:"
         write(6,*) "Invalid index, beta =",beta
         call exit(-10)
      endif
*
*     interpolate
*
      n = inter_degree(0)
      do alpha=0,beta
         if(xg(0,alpha).gt.x) goto 11
      enddo
 11   lx = alpha - 1
      ux = lx + n + 1
      ExternalEvolutionOperator = 0d0
      if(Bs2Bs.eq."Ev2Ev")then
         do alpha=lx,ux
            wg = w_int_gen(n,alpha,x)
            ExternalEvolutionOperator = ExternalEvolutionOperator 
     1                                + wg * Ev2EvQCD(0,i,j,alpha,beta)
         enddo
      elseif(Bs2Bs.eq."Ev2Ph")then
         do alpha=lx,ux
            wg = w_int_gen(n,alpha,x)
            ExternalEvolutionOperator = ExternalEvolutionOperator 
     1                                + wg * Ev2PhQCD(0,i,j,alpha,beta)
         enddo
      elseif(Bs2Bs.eq."Ph2Ph")then
         do alpha=lx,ux
            wg = w_int_gen(n,alpha,x)
            ExternalEvolutionOperator = ExternalEvolutionOperator 
     1                                + wg * Ph2PhQCD(0,i,j,alpha,beta)
         enddo
      endif
*
      return
      end
