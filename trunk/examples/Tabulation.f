************************************************************************
*
*     Tabulation.f:
*
*     Example program used for the LH benchmark.
*
************************************************************************
      program Tabulation
*
      implicit none
*
      integer ilha
      double precision Q0,Q
      double precision Q02,Q2
      double precision AlphaQCD,AlphaQED
      double precision xPDF,xgamma
      double precision xPDFj,xgammaj
      double precision eps
      double precision xlha(11)

      parameter(eps=1d-10)
      data xlha / 1d-7, 1d-6, 1d-5, 1d-4, 1d-3, 1d-2,
     1            1d-1, 3d-1, 5d-1, 7d-1, 9d-1 /

      integer n,alpha,beta,i,j,k
      double precision xext(0:200)!,step
      double precision M(0:14*14*(201)*(201)-1)
      double precision f(-6:6,0:200),f0(-6:6,0:200),xpd(-6:6),fin(-6:6)
      double precision w_int_herm
*
*     Some examples ...
*
c      call LockGrids(.true.)
c      call EnableEvolutionOperator(.true.)
c      call SetFFNS(3)
c      call SetTheory("QECDS")
      call SetPerturbativeOrder(0)
c      call SetPDFSet("NNPDF23_nlo_as_0119_qed.LHgrid")
c      call SetPDFSet("MRST2004qed.LHgrid")
c      call SetNumberOfGrids(1)
c      call SetGridParameters(1,30,3,1d-5)
c      call SetGridParameters(2,30,3,2d-1)
c      call SetGridParameters(3,30,3,8d-1)
*
*     Initializes integrals on the grids
*
      call InitializeAPFEL
*
*     Evolve PDFs on the grids
*
      write(6,*) "Enter initial and final scale in GeV^2"
      read(5,*) Q02,Q2
*
      Q0 = dsqrt(Q02) - eps
      Q  = dsqrt(Q2)
      call EvolveAPFEL(Q0,Q)
*
*     Tabulate PDFs for the LHA x values
*
      write(6,*) "alpha_QCD(mu2F) =",AlphaQCD(Q)
      write(6,*) "alpha_QED(mu2F) =",AlphaQED(Q)
      write(6,*) "  "
*
      write(6,*) "Standard evolution:"
      write(6,'(a5,2a12,a14,a10,2a12)') "x",
     1         "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon","photon"
      do ilha=3,11
         write(6,'(es7.1,6es12.4)') 
     1         xlha(ilha),
     2         xPDF(2,xlha(ilha)) - xPDF(-2,xlha(ilha)),
     3         xPDF(1,xlha(ilha)) - xPDF(-1,xlha(ilha)),
     4         2d0 * ( xPDF(-1,xlha(ilha)) + xPDF(-2,xlha(ilha)) ),
     5         xPDF(4,xlha(ilha)) + xPDF(-4,xlha(ilha)),
     6         xPDF(0,xlha(ilha)),
     7         xgamma(xlha(ilha))
      enddo
      write(*,*) "  "
*
      write(6,*) "Evolution from the joint grid:"
      write(6,'(a5,2a12,a14,a10,2a12)') "x",
     1         "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon","photon"
      do ilha=3,11
         write(6,'(es7.1,6es12.4)') 
     1         xlha(ilha),
     2         xPDFj(2,xlha(ilha)) - xPDFj(-2,xlha(ilha)),
     3         xPDFj(1,xlha(ilha)) - xPDFj(-1,xlha(ilha)),
     4         2d0 * ( xPDFj(-1,xlha(ilha)) + xPDFj(-2,xlha(ilha)) ),
     5         xPDFj(4,xlha(ilha)) + xPDFj(-4,xlha(ilha)),
     6         xPDFj(0,xlha(ilha)),
     7         xgammaj(xlha(ilha))
      enddo
      write(*,*) "  "
*
*     Generate (or read) an external grid
*
c      n = 40
c      xext(0) = 1d-5
c*     log
c      step = dlog(1d0/xext(0)) / dble(n) 
c      do alpha=1,n
c         xext(alpha) = xext(alpha-1) * dexp(step)
c      enddo
c*     lin
c      step = ( 1d0 - xext(0) ) / dble(n) 
c      do alpha=1,n
c         xext(alpha) = xext(alpha-1) + step
c      enddo
*     External file
      open(unit=12,file="APPLgrid.dat",status="unknown")
      read(12,*) n
      do alpha=0,n
         read(12,*) xext(alpha)
      enddo
      close(12)
*
*     Initial PDFs
*
      do alpha=0,n
         call toyLHPDFs(xext(alpha),xpd)
         do i=-6,6
            f0(i,alpha) = xpd(i)
         enddo
      enddo
*
      call ExternalEvolutionOperator(Q0,Q,n,xext,M)
*
*     Evolve PDFs with "M"
*
      do alpha=0,n
         do i=-6,6
            f(i,alpha) = 0d0
            do beta=0,n
               do j=-6,6
                  k = ( 7 + i ) + 14 * ( ( 7 + j ) 
     1              + 14 * ( alpha + ( n + 1 ) * beta ) )
                  f(i,alpha) = f(i,alpha) + M(k) * f0(j,beta)
               enddo
            enddo
         enddo
      enddo
*
      write(6,*) "Evolution with 'ExternalEvolutionOperator':"
      write(6,'(a5,2a12,a14,a10,2a12)') "x",
     1         "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon","photon"
      do ilha=3,11
         do i=-6,6
            fin(i) = 0d0
            do alpha=0,n
               fin(i) = fin(i) 
     1                + w_int_herm(n,xext,alpha,xlha(ilha)) * f(i,alpha)
            enddo
            fin(i) = 100d0 * ( fin(i) - xPDF(i,xlha(ilha)) ) 
     1             / xPDF(i,xlha(ilha))
         enddo
         write(6,'(es7.1,6es12.4)') 
     1         xlha(ilha),
     2         fin(2) - fin(-2),
     3         fin(1) - fin(-1),
     4         2d0 * ( fin(-1) + fin(-2) ),
     5         fin(4) + fin(-4),
     6         fin(0),
     7         0d0         
      enddo
      write(*,*) "  "
*
      end
