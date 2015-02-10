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
      function fxpdfall(x,xf)
      double precision x, xf(-6:6), fxpdfall
      fxpdfall = xPDFall(x,xf)
      return
      end

ccccccccccccc
      function fexternalevolutionoperator(fname,i,j,x,beta)
      integer i,j,beta
      double precision x,fexternalevolutionoperator
      character fname*(*)
      fexternalevolutionoperator = 
     1     ExternalEvolutionOperator(fname,i,j,alpha,x)
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
      subroutine fsetfastevolution(fe)
      logical fe
      call SetFastEvolution(fe)
      end subroutine fsetfastevolution

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

cccc Functions for the new DIS module ccccc
      subroutine finitializeapfel_dis
      call InitializeAPFEL_DIS
      end subroutine finitializeapfel_dis

ccccccccccccc
      subroutine fcomputestructurefunctionsapfel(Q0,Q)
      double precision Q,Q0
      call ComputeStructureFunctionsAPFEL(Q0,Q)
      end subroutine fcomputestructurefunctionsapfel

ccccccccccccc
      subroutine fsetmassscheme(ms)
      character ms*(*)
      call SetMassScheme(ms)
      end subroutine fsetmassscheme

ccccccccccccc
      subroutine fsetpolarizationdis(pol)
      double precision pol
      call SetPolarizationDIS(pol)
      end subroutine fsetpolarizationdis

ccccccccccccc
      subroutine fsetprocessdis(pr)
      character pr*(*)
      call SetProcessDIS(pr)
      end subroutine fsetprocessdis

ccccccccccccc
      subroutine fsetprojectiledis(lept)
      character lept*(*)
      call SetProjectileDIS(lept)
      end subroutine fsetprojectiledis

ccccccccccccc
      subroutine fsettargetdis(tar)
      character tar*(*)
      call SetTargetDIS(tar)
      end subroutine fsettargetdis

ccccccccccccc
      function fexternaldisoperator(SF,ihq,i,x,beta)
      integer ihq,i
      double precision x
      integer beta
      character SF*(*)
      double precision fexternaldisoperator
      fexternaldisoperator = ExternalDISOperator(SF,ihq,i,x,beta)
      return
      end

ccccccccccccc
      function ff2light(x)
      double precision x,ff2light
      ff2light = F2light(x)
      return
      end

ccccccccccccc
      function ff2charm(x)
      double precision x,ff2charm
      ff2charm = F2charm(x)
      return
      end

ccccccccccccc
      function ff2bottom(x)
      double precision x,ff2bottom
      ff2bottom = F2bottom(x)
      return
      end

ccccccccccccc
      function ff2top(x)
      double precision x,ff2top
      ff2top = F2top(x)
      return
      end

ccccccccccccc
      function ff2total(x)
      double precision x,ff2total
      ff2total = F2total(x)
      return
      end

ccccccccccccc
      function ffllight(x)
      double precision x,ffllight
      ffllight = FLlight(x)
      return
      end

ccccccccccccc
      function fflcharm(x)
      double precision x,fflcharm
      fflcharm = FLcharm(x)
      return
      end

ccccccccccccc
      function fflbottom(x)
      double precision x,fflbottom
      fflbottom = FLbottom(x)
      return
      end

ccccccccccccc
      function ffltop(x)
      double precision x,ffltop
      ffltop = FLtop(x)
      return
      end

ccccccccccccc
      function ffltotal(x)
      double precision x,ffltotal
      ffltotal = FLtotal(x)
      return
      end

ccccccccccccc
      function ff3light(x)
      double precision x,ff3light
      ff3light = F3light(x)
      return
      end

ccccccccccccc
      function ff3charm(x)
      double precision x,ff3charm
      ff3charm = F3charm(x)
      return
      end

ccccccccccccc
      function ff3bottom(x)
      double precision x,ff3bottom
      ff3bottom = F3bottom(x)
      return
      end

ccccccccccccc
      function ff3top(x)
      double precision x,ff3top
      ff3top = F3top(x)
      return
      end

ccccccccccccc
      function ff3total(x)
      double precision x,ff3total
      ff3total = F3total(x)
      return
      end

ccccccccccccc
      subroutine fsetzmass(massz)
      double precision massz
      call SetZMass(massz)
      end subroutine fsetzmass

ccccccccccccc
      subroutine fsetwmass(massw)
      double precision massw
      call SetWMass(massw)
      end subroutine fsetwmass

ccccccccccccc
      subroutine fsetprotonmass(massp)
      double precision massp
      call SetProtonMass(massp)
      end subroutine fsetprotonmass

ccccccccccccc
      subroutine fsetsinthetaw(sw)
      double precision sw
      call SetSinThetaW(sw)
      end subroutine fsetsinthetaw

ccccccccccccc
      subroutine fsetckm(vud,vus,vub,vcd,vcs,vcb,vtd,vts,vtb)
      double precision vud,vus,vub,vcd,vcs,vcb,vtd,vts,vtb
      call SetCKM(vud,vus,vub,vcd,vcs,vcb,vtd,vts,vtb)
      end subroutine fsetckm

ccccccccccccc
      subroutine fsetgfermi(gf)
      double precision gf
      call SetGFermi(gf)
      end subroutine fsetgfermi

ccccccccccccc
      function ffksimulator(obs,x,q,y,i,beta)
      integer i,beta
      double precision x,q,y
      character obs*(*)
      double precision ffksimulator
      ffksimulator = FKSimulator(obs,x,q,y,i,beta)
      return
      end

ccccccccccccc
      subroutine fsetfkobservable(obs)
      character obs*(*)
      call SetFKObservable(obs)
      return
      end
