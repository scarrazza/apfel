C     The APFEL Fortran wrapper

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
      subroutine fselectcharge(selch)
      character selch*(*)
      call SelectCharge(selch)
      end subroutine fselectcharge

ccccccccccccc      
      subroutine fsetrenqratio(ratioR)
      double precision ratioR
      call SetRenQRatio(ratioR)
      end subroutine fsetrenqratio

ccccccccccccc      
      subroutine fsetfacqratio(ratioF)
      double precision ratioF
      call SetFacQRatio(ratioF)
      end subroutine fsetfacqratio

ccccccccccccc      
      subroutine fenabledynamicalscalevariations(dsv)
      logical dsv
      call EnableDynamicalScaleVariations(dsv)
      end subroutine fenabledynamicalscalevariations

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
      subroutine fsetsin2thetaw(sw)
      double precision sw
      call SetSin2ThetaW(sw)
      end subroutine fsetsin2thetaw

ccccccccccccc
      subroutine fsetckm(vud,vus,vub,vcd,vcs,vcb,vtd,vts,vtb)
      double precision vud,vus,vub,vcd,vcs,vcb,vtd,vts,vtb
      call SetCKM(vud,vus,vub,vcd,vcs,vcb,vtd,vts,vtb)
      end subroutine fsetckm

ccccccccccccc
      subroutine fsetpropagatorcorrection(dr)
      double precision dr
      call SetPropagatorCorrection(dr)
      end subroutine fsetpropagatorcorrection

ccccccccccccc
      subroutine fsetewcouplings(vd,vu,ad,au)
      double precision vd,vu,ad,au
      call SetEWCouplings(vd,vu,ad,au)
      end subroutine fsetewcouplings

ccccccccccccc
      subroutine fsetgfermi(gf)
      double precision gf
      call SetGFermi(gf)
      end subroutine fsetgfermi

ccccccccccccc
      function fgetzmass()
      double precision fgetzmass
      fgetzmass = GetZMass()
      return
      end

ccccccccccccc
      function fgetwmass()
      double precision fgetwmass
      fgetwmass = GetWMass()
      return
      end

ccccccccccccc
      function fgetprotonmass()
      double precision fgetprotonmass
      fgetprotonmass = GetProtonMass()
      return
      end

ccccccccccccc
      function fgetsin2thetaw()
      double precision fgetsin2thetaw
      fgetsin2thetaw = GetSin2ThetaW()
      return
      end

ccccccccccccc
      function fgetckm(u,d)
      integer u,d
      double precision fgetckm
      fgetckm = GetCKM(u,d)
      return
      end

ccccccccccccc
      function fgetgfermi()
      double precision fgetgfermi
      fgetgfermi = GetGFermi()
      return
      end

ccccccccccccc
      function fgetsiatotalcrosssection(pto,Q)
      integer pto
      double precision Q
      double precision fgetsiatotalcrosssection
      fgetsiatotalcrosssection = GetSIATotalCrossSection(pto,Q)
      return
      end

ccccccccccccc
      subroutine fenabletargetmasscorrections(tc)
      logical tc
      call EnableTargetMassCorrections(tc)
      end subroutine fenabletargetmasscorrections

ccccccccccccc
      subroutine fenabledampingfonll(df)
      logical df
      call EnableDampingfonll(df)
      end subroutine fenabledampingfonll

ccccccccccccc
      function ffksimulator(x,q,y,i,beta)
      integer i,beta
      double precision x,q,y
      double precision ffksimulator
      ffksimulator = FKSimulator(x,q,y,i,beta)
      return
      end

ccccccccccccc
      subroutine fsetfkobservable(obs)
      character obs*(*)
      call SetFKObservable(obs)
      return
      end

ccccccccccccc
      subroutine fgetfkobservable
      call GetFKObservable
      end subroutine fgetfkobservable

ccccccccccccc
      function ffkobservables(x,q,y)
      double precision x,q,y
      double precision ffkobservables
      ffkobservables = FKObservables(x,q,y)
      return
      end

ccccccccccccc
      subroutine fcomputefktables(inputfile,outputpath,Q0,flmap)
      double precision Q0
      character inputfile*(*)
      character outputpath*(*)
      integer flmap(0:13,0:13)
      call ComputeFKTables(inputfile,outputpath,Q0,flmap)
      end subroutine fcomputefktables

ccccccccccccc
      subroutine fcomputehardcrosssectionsdy(inputfile,outputfile)      
      character inputfile*(*)
      character outputfile*(*)
      call ComputeHardCrossSectionsDY(inputfile,outputfile)
      end subroutine fcomputehardcrosssectionsdy
