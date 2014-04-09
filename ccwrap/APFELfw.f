!     The APFEL Fortran wrapper

ccccccccccccc
      subroutine finitializeapfel
      call InitializeAPFEL
      end subroutine finitializeapfel

ccccccccccccc      
      subroutine fevolveapfel(Q0,Q)
      double precision Q0, Q
      call EvolveAPFEL(Q0,Q)
      end subroutine

ccccccccccccc      
      function fxpdf(i,x)
      integer i
      double precision x, fxpdf
      fxpdf = xPDF(i,x)
      return
      end

ccccccccccccc      
      function fxpdfj(i,x)
      integer i
      double precision x, fxpdfj
      fxpdfj = xPDFj(i,x)
      return
      end

ccccccccccccc      
      function fdxpdf(i,x)
      integer i
      double precision x, fdxpdf
      fdxpdf = dxPDF(i,x)
      return
      end

ccccccccccccc      
      function fxgamma(x)
      double precision x, fxgamma
      fxgamma = xgamma(x)
      return
      end

ccccccccccccc      
      function fxgammaj(x)
      double precision x, fxgammaj
      fxgammaj = xgammaj(x)
      return
      end

ccccccccccccc      
      function fdxgamma(x)
      double precision x, fdxgamma
      fdxgamma = dxgamma(x)
      return
      end

ccccccccccccc      
      subroutine flhapdfgrid(Nrep,Qin,fname)
      integer Nrep
      double precision Qin
      character fname*(*)
      call LHAPDFgrid(Nrep, Qin, fname)
      end subroutine flhapdfgrid

ccccccccccccc      
      subroutine flhapdfgridderivative(Nrep,fname)
      integer Nrep
      character fname*(*)
      call LHAPDFgridDerivative(Nrep, fname)
      end subroutine flhapdfgridderivative

ccccccccccccc
      subroutine fexternalevolutionoperator(q0,q,n,xext,m)
      integer n
      double precision q0,q
      double precision xext(0:n)
      double precision m(0:14*14*(n+1)*(n+1)-1)
      call ExternalEvolutionOperator(q0,q,n,xext,m)
      end subroutine fexternalevolutionoperator

ccccccccccccc      
      function falphaqcd(Q)
      double precision Q, falphaqcd
      falphaqcd = AlphaQCD(Q)
      return
      end

ccccccccccccc      
      function falphaqed(Q)
      double precision Q, falphaqed
      falphaqed = AlphaQED(Q)
      return
      end

ccccccccccccc      
      function fnpdf(i,N)
      integer i,N
      double precision fnpdf
      fnpdf = NPDF(i,N)
      return
      end

ccccccccccccc      
      function fngamma(N)
      integer N
      double precision fngamma
      fngamma = Ngamma(N)
      return
      end

ccccccccccccc      
      function flumi(i,j,S)
      integer i,j
      double precision S,flumi
      flumi = LUMI(i,j,S)
      return
      end

ccccccccccccc      
      subroutine fcleanup
      call CleanUp
      end subroutine fcleanup

ccccccccccccc      
      subroutine fenablewelcomemessage(wc)
      logical wc
      call EnableWelcomeMessage(wc)
      end subroutine fenablewelcomemessage

ccccccccccccc      
      subroutine fenableevolutionoperator(eo)
      logical eo
      call EnableEvolutionOperator(eo)
      end subroutine fenableevolutionoperator

ccccccccccccc      
      subroutine flockgrids(lg)
      logical lg
      call LockGrids(lg)
      end subroutine flockgrids

ccccccccccccc      
      subroutine flsettimelikeevolution(tl)
      logical tl
      call SetTimeLikeEvolution(tl)
      end subroutine flsettimelikeevolution

ccccccccccccc      
      function fheavyquarkmass(i,Q)
      integer i
      double precision Q,fheavyquarkmass
      fheavyquarkmass = HeavyQuarkMass(i,Q)
      return
      end

ccccccccccccc      
      subroutine fsetalphaqcdref(alpharef,Qref)
      double precision alpharef,Qref      
      call SetAlphaQCDRef(alpharef,Qref)
      end subroutine fsetalphaqcdref

ccccccccccccc      
      subroutine fsetalphaqedref(alpharef,Qref)
      double precision alpharef,Qref      
      call SetAlphaQEDRef(alpharef,Qref)
      end subroutine fsetalphaqedref

ccccccccccccc      
      subroutine fsetalphaevolution(evol)
      character evol*(*)
      call SetAlphaEvolution(evol)
      end subroutine fsetalphaevolution

ccccccccccccc      
      subroutine fsetlambdaqcdref(lambdaref,nref)
      integer nref
      double precision lambdaref
      call SetLambdaQCDRef(lambdaref,nref)
      end subroutine fsetlambdaqcdref

ccccccccccccc      
      subroutine fsetqlimits(Qmin,Qmax)
      double precision Qmin,Qmax 
      call SetQLimits(Qmin,Qmax)
      end subroutine fsetqlimits

ccccccccccccc      
      subroutine fsetffns(nfl)
      integer nfl
      call SetFFNS(nfl)
      end subroutine fsetffns

ccccccccccccc      
      subroutine fsetgridparameters(i,np,deg,x)
      integer i, np, deg
      double precision x
      call SetGridParameters(i,np,deg,x)
      end subroutine fsetgridparameters

ccccccccccccc      
      subroutine fsetmaxflavouralpha(nf)
      integer nf
      call SetMaxFlavourAlpha(nf)
      end subroutine fsetmaxflavouralpha

ccccccccccccc      
      subroutine fsetmaxflavourpdfs(nf)
      integer nf
      call SetMaxFlavourPDFs(nf)
      end subroutine fsetmaxflavourpdfs

ccccccccccccc      
      subroutine fsetmsbarmasses(mc,mb,mt)
      double precision mc, mb, mt
      call SetMSbarMasses(mc,mb,mt)
      end subroutine fsetmsbarmasses

ccccccccccccc      
      subroutine fsetnumberofgrids(n)
      integer n
      call SetNumberOfGrids(n)
      end subroutine fsetnumberofgrids

ccccccccccccc      
      subroutine fsetpdfset(name)
      character name*(*)
      call SetPDFSet(name)
      end subroutine fsetpdfset

ccccccccccccc      
      subroutine fsetperturbativeorder(pto)
      integer pto
      call SetPerturbativeOrder(pto)
      end subroutine fsetperturbativeorder

ccccccccccccc      
      subroutine fsetpolemasses(mc,mb,mt)
      double precision mc, mb, mt
      call SetPoleMasses(mc,mb,mt)
      end subroutine fsetpolemasses

ccccccccccccc      
      subroutine fsetrenfacratio(ratio)
      double precision ratio
      call SetRenFacRatio(ratio)
      end subroutine fsetrenfacratio

ccccccccccccc      
      subroutine fsetreplica(nr)
      integer nr
      call SetReplica(nr)
      end subroutine fsetreplica

ccccccccccccc      
      subroutine fsettheory(theory)
      character theory*(*)
      call SetTheory(theory)
      end subroutine fsettheory

ccccccccccccc      
      subroutine fsetvfns 
      call SetVFNS
      end subroutine fsetvfns

ccccccccccccc
      subroutine fdisxsec(x,qi,qf,y,pol,proc,scheme,pto,pdfset,irep,
     1                    target,proj,F2,F3,FL,SIGMA)
      double precision x,qi,qf,y,pol
      character  proc*(*)
      character  scheme*(*)
      character pdfset*(*)
      character  target*(*)
      character proj*(*)
      integer pto,irep
      double precision F2(3:7),F3(3:7),FL(3:7),SIGMA(3:7)
      call DIS_xsec(x,qi,qf,y,pol,proc,scheme,pto,pdfset,irep,
     1              target,proj,F2,F3,FL,SIGMA)
      end subroutine fdisxsec

