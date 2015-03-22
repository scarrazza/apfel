************************************************************************
*
*     TabulationStep.f:
*
*     Example program used for the LH benchmark but using the stepping
*     of the evolution.
*
************************************************************************
      program TabulationStep
*
      implicit none
*
      integer ilha
      integer i,n
      integer ipt
      double precision Q0,Q,Qi,Qf
      double precision Q02,Q2
      double precision xPDF,xgamma
      double precision eps
      double precision delta
      double precision xlha(11)

      parameter(eps=1d-10)
      data xlha / 1d-7, 1d-6, 1d-5, 1d-4, 1d-3, 1d-2,
     1            1d-1, 3d-1, 5d-1, 7d-1, 9d-1 /
*
*     Evolve PDFs on the grids
*
      write(6,*) "Enter initial and final scale in GeV^2"
      read(5,*) Q02,Q2
*
      Q0 = dsqrt(Q02) - eps
      Q  = dsqrt(Q2)
*
      ipt = 1
*
*     QECD
*
      call SetTheory("QECDP")
      call SetFastEvolution(.false.)
      call SetQLimits(0.5d0,20000d0)
      call SetPerturbativeOrder(ipt)
      call InitializeAPFEL
      do n=1,100,99
         write(6,*) "QECDP with",n," steps"
         write(6,*) "  "
         call SetPDFSet("ToyLH")
         delta = dlog(Q/Q0) / dble(n)     
         Qi = Q0
         do i=1,n
            Qf = Qi * dexp(delta)
            call EvolveAPFEL(Qi,Qf)
            call SetPDFSet("apfel")
            Qi = Qf
         enddo
         write(6,'(a5,2a12,a14,a10,2a12)') "x",
     1        "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon","photon"
         do ilha=3,11
            write(6,'(es7.1,6es12.4)') 
     1           xlha(ilha),
     2           xPDF(2,xlha(ilha)) - xPDF(-2,xlha(ilha)),
     3           xPDF(1,xlha(ilha)) - xPDF(-1,xlha(ilha)),
     4           2d0 * ( xPDF(-1,xlha(ilha)) + xPDF(-2,xlha(ilha)) ),
     5           xPDF(4,xlha(ilha)) + xPDF(-4,xlha(ilha)),
     6           xPDF(0,xlha(ilha)),
     7           xgamma(xlha(ilha))
         enddo
         write(*,*) "  "
      enddo
      call CleanUp
*
*     QCED
*
      call SetTheory("QCEDP")
      call EnableWelcomeMessage(.false.)
      call SetFastEvolution(.false.)
      call SetQLimits(0.5d0,20000d0)
      call SetPerturbativeOrder(ipt)
      call InitializeAPFEL
      do n=1,100,99
         write(6,*) "QCEDP with",n," steps"
         write(6,*) "  "
         call SetPDFSet("ToyLH")
         delta = dlog(Q/Q0) / dble(n)     
         Qi = Q0
         do i=1,n
            Qf = Qi * dexp(delta)
            call EvolveAPFEL(Qi,Qf)
            call SetPDFSet("apfel")
            Qi = Qf
         enddo
         write(6,'(a5,2a12,a14,a10,2a12)') "x",
     1        "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon","photon"
         do ilha=3,11
            write(6,'(es7.1,6es12.4)') 
     1           xlha(ilha),
     2           xPDF(2,xlha(ilha)) - xPDF(-2,xlha(ilha)),
     3           xPDF(1,xlha(ilha)) - xPDF(-1,xlha(ilha)),
     4           2d0 * ( xPDF(-1,xlha(ilha)) + xPDF(-2,xlha(ilha)) ),
     5           xPDF(4,xlha(ilha)) + xPDF(-4,xlha(ilha)),
     6           xPDF(0,xlha(ilha)),
     7           xgamma(xlha(ilha))
         enddo
         write(*,*) "  "
      enddo
      call CleanUp
*
*     QavD
*
      call SetTheory("QavDP")
      call EnableWelcomeMessage(.false.)
      call SetFastEvolution(.false.)
      call SetQLimits(0.5d0,20000d0)
      call SetPerturbativeOrder(ipt)
      call InitializeAPFEL
      do n=1,100,99
         write(6,*) "QavDP with",n," steps"
         write(6,*) "  "
         call SetPDFSet("ToyLH")
         delta = dlog(Q/Q0) / dble(n)     
         Qi = Q0
         do i=1,n
            Qf = Qi * dexp(delta)
            call EvolveAPFEL(Qi,Qf)
            call SetPDFSet("apfel")
            Qi = Qf
         enddo
         write(6,'(a5,2a12,a14,a10,2a12)') "x",
     1        "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon","photon"
         do ilha=3,11
            write(6,'(es7.1,6es12.4)') 
     1           xlha(ilha),
     2           xPDF(2,xlha(ilha)) - xPDF(-2,xlha(ilha)),
     3           xPDF(1,xlha(ilha)) - xPDF(-1,xlha(ilha)),
     4           2d0 * ( xPDF(-1,xlha(ilha)) + xPDF(-2,xlha(ilha)) ),
     5           xPDF(4,xlha(ilha)) + xPDF(-4,xlha(ilha)),
     6           xPDF(0,xlha(ilha)),
     7           xgamma(xlha(ilha))
         enddo
         write(*,*) "  "
      enddo
      call CleanUp
*
*     QUniD
*
      call SetTheory("QUniD")
      call EnableWelcomeMessage(.false.)
c      call SetFastEvolution(.false.)
      call SetQLimits(0.5d0,20000d0)
      call SetPerturbativeOrder(ipt)
      call InitializeAPFEL
      do n=1,100,99
         write(6,*) "QUniD with",n," steps"
         write(6,*) "  "
         call SetPDFSet("ToyLH")
         delta = dlog(Q/Q0) / dble(n)     
         Qi = Q0
         do i=1,n
            Qf = Qi * dexp(delta)
            call EvolveAPFEL(Qi,Qf)
            call SetPDFSet("apfel")
            Qi = Qf
         enddo
         write(6,'(a5,2a12,a14,a10,2a12)') "x",
     1        "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon","photon"
         do ilha=3,11
            write(6,'(es7.1,6es12.4)') 
     1           xlha(ilha),
     2           xPDF(2,xlha(ilha)) - xPDF(-2,xlha(ilha)),
     3           xPDF(1,xlha(ilha)) - xPDF(-1,xlha(ilha)),
     4           2d0 * ( xPDF(-1,xlha(ilha)) + xPDF(-2,xlha(ilha)) ),
     5           xPDF(4,xlha(ilha)) + xPDF(-4,xlha(ilha)),
     6           xPDF(0,xlha(ilha)),
     7           xgamma(xlha(ilha))
         enddo
         write(*,*) "  "
      enddo
      call CleanUp
*
      end
