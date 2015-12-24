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
      subroutine fderiveapfel(Q)
      double precision Q
      call DeriveAPFEL(Q)
      end subroutine

ccccccccccccc      
      subroutine fcachepdfsapfel(Q0)
      double precision Q0
      call CachePDFsAPFEL(Q0)
      end subroutine

ccccccccccccc      
      function fxpdf(i,x)
      integer i
      double precision x, fxpdf
      fxpdf = xPDF(i,x)
      return
      end

ccccccccccccc      
      function fxpdfxq(i,x,Q)
      integer i
      double precision x, Q, fxpdfxq
      fxpdfxq = xPDFxQ(i,x,Q)
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
      function fxpdfall(x,xf)
      double precision x, xf(-6:6), fxpdfall
      fxpdfall = xPDFall(x,xf)
      return
      end

ccccccccccccc      
      function fxlepton(i,x)
      integer i
      double precision x, fxlepton
      fxlepton = xLepton(i,x)
      return
      end

ccccccccccccc      
      function fxleptonj(i,x)
      integer i
      double precision x, fxleptonj
      fxleptonj = xLeptonj(i,x)
      return
      end

ccccccccccccc
      function fexternalevolutionoperator(fname,i,j,x,beta)
      integer i,j,beta
      double precision x,fexternalevolutionoperator
      double precision ExternalEvolutionOperator
      character fname*(*)
      fexternalevolutionoperator = 
     1     ExternalEvolutionOperator(fname,i,j,x,beta)
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
      function fxgrid(alpha)
      integer alpha
      double precision fxgrid
      fxgrid = xGrid(alpha)
      return
      end

ccccccccccccc      
      function fnintervals()
      integer nintervals
      fnintervals = nIntervals()
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
      subroutine fenableleptonevolution(le)
      logical le
      call EnableLeptonEvolution(le)
      end subroutine fenableleptonevolution

ccccccccccccc      
      subroutine flockgrids(lg)
      logical lg
      call LockGrids(lg)
      end subroutine flockgrids

ccccccccccccc      
      subroutine fsettimelikeevolution(tl)
      logical tl
      call SetTimeLikeEvolution(tl)
      end subroutine fsettimelikeevolution

ccccccccccccc      
      subroutine fsetpolarizedevolution(polev)
      logical polev
      call SetPolarizedEvolution(polev)
      end subroutine fsetpolarizedevolution

ccccccccccccc      
      subroutine fsetfastevolution(fe)
      logical fe
      call SetFastEvolution(fe)
      end subroutine fsetfastevolution

ccccccccccccc      
      subroutine fenablemassrunning(mr)
      logical mr
      call EnableMassRunning(mr)
      end subroutine fenablemassrunning

ccccccccccccc      
      subroutine fsetsmallxresummation(sx,la)
      logical sx
      character la*(*)
      call SetSmallxResummation(sx,la)
      end subroutine fsetsmallxresummation

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
      subroutine fsetpdfevolution(evolp)
      character evolp*(*)
      call SetPDFEvolution(evolp)
      end subroutine fsetpdfevolution

ccccccccccccc      
      subroutine fsetepsilontruncation(eps)
      double precision eps
      call SetEpsilonTruncation(eps)
      end subroutine fsetepsilontruncation

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
      subroutine fsetqgridparameters(npq,degq)
      integer npq, degq
      call SetQGridParameters(npq,degq)
      end subroutine fsetqgridparameters

ccccccccccccc
      subroutine fsetexternalgrid(i,np,deg,x)
      integer i, np, deg
      double precision x(0:np)
      call SetExternalGrid(i,np,deg,x)
      end subroutine fsetexternalgrid

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
      subroutine fsetmassscalereference(Qc,Qb,Qt)
      double precision Qc, Qb, Qt
      call SetMassScaleReference(Qc,Qb,Qt)
      end subroutine fsetmassscalereference

ccccccccccccc      
      subroutine fsetmassmatchingscales(kmc,kmb,kmt)
      double precision kmc, kmb, kmt
      call SetMassMatchingScales(kmc,kmb,kmt)
      end subroutine fsetmassmatchingscales

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
      function fgetperturbativeorder()
      integer fgetperturbativeorder
      integer GetPerturbativeOrder
      fgetperturbativeorder = GetPerturbativeOrder()
      return
      end

ccccccccccccc      
      function fgetmuf()
      double precision fgetmuf
      double precision GetMuF
      fgetmuf = GetMuF()
      return
      end

ccccccccccccc      
      function fgetmuf0()
      double precision fgetmuf0
      double precision GetMuF0
      fgetmuf0 = GetMuF0()
      return
      end

ccccccccccccc      
      subroutine fsetpolemasses(mc,mb,mt)
      double precision mc, mb, mt
      call SetPoleMasses(mc,mb,mt)
      end subroutine fsetpolemasses

ccccccccccccc
      subroutine fsettaumass(masst)
      double precision masst
      call SetTauMass(masst)
      end subroutine fsettaumass

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
      subroutine flistfunctions
      call ListFunctions
      end subroutine flistfunctions

ccccccccccccc
      function fcheckapfel()
      logical fcheckapfel, CheckAPFEL
      fcheckapfel = CheckAPFEL()
      return
      end
