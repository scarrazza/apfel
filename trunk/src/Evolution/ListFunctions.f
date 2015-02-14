************************************************************************
*
*     ListFunctions.f:
*
*     This subroutine lists all the functions available in APFEL.
*     It needs to be update constantly!
*
************************************************************************
      subroutine ListFunctions
*
      implicit none
*
      call WelcomeMessage
*
      write(6,*) "List of the functions available in APFEL:"
      write(6,*) "  "
c      write(6,*) "Initialize the library:"
      write(6,*) "- InitializeAPFEL()"
c      write(6,*) "Compute evolution:"
      write(6,*) "- EvolveAPFEL(double Q0, double Q)"
c      write(6,*) "Computes the derivative of PDFs:"
      write(6,*) "- DeriveAPFEL(double Q)"
c      write(6,*) "Returns x*PDF
      write(6,*) "- xPDF(int i, double x)"
c      write(6,*) "Returns the derivative of x*PDF
      write(6,*) "- dxPDF(int i, double x)"
c      write(6,*) "Returns x*PDF on the joint grid
      write(6,*) "- xPDFj(int i, double x)"
c      write(6,*) "Returns x*gamma
      write(6,*) "- xgamma(double x)"
c      write(6,*) "Returns x*gamma on the joint grid
      write(6,*) "- xgammaj(double x)"
c      write(6,*) "Returns the derivative of x*gamma
      write(6,*) "- dxgamma(double x)"
c      write(6,*) "Returns all PDFs on the joint grid
      write(6,*) "- xPDFall(double x, double *xf)"
c      write(6,*) "External Evolution Operator
      write(6,*) "- ExternalEvolutionOperator(string fname,",
     1           " int i, int j, double x, int beta)"
c      write(6,*) "Builds the *.LHgrid output file
      write(6,*) "- LHAPDFgrid(int Nrep, double Qin, string",
     1           " fname)"
c      write(6,*) "Builds the *.LHgrid output file with the derivative of the input set
      write(6,*) "- LHAPDFgridDerivative(int Nrep, string",
     1           " fname)"
c      write(6,*) "Returns the value of alpha_QCD at the given scale
      write(6,*) "- AlphaQCD(double Q)"
c      write(6,*) "Returns the value of alpha_QED at the given scale
      write(6,*) "- AlphaQED(double Q)"
c      write(6,*) "Returns Heavy Quark Masses
      write(6,*) "- HeavyQuarkMass(int,double)"
c      write(6,*) "Returns the N-th Mellin moment of the i-th PDF 
c      write(6,*) "in the physical basis at the final scale
      write(6,*) "- NPDF(int i, int N)"
c      write(6,*) "Returns the N-th Mellin moment of the photon PDF 
c      write(6,*) "in the physical basis at the final scale
      write(6,*) "- Ngamma(int N)"
c      write(6,*) "Returns the Luminosity 
      write(6,*) "- LUMI(int i, int j, double S)"
c      write(6,*) "Returns the joint x-space grid 
      write(6,*) "- xGrid(int alpha)"
c      write(6,*) "Returns the number of intervals of the joint x-space grid
      write(6,*) "- nIntervals()"
c      write(6,*) "Gets APFEL version
      write(6,*) "- GetVersion()"
c      write(6,*) "- /*
c      write(6,*) "- * Set Methods
c      write(6,*) "- */
c      write(6,*) "Cleans up parameters
      write(6,*) "- CleanUp()"
c      write(6,*) "Enables welcome message
      write(6,*) "- EnableWelcomeMessage(bool)"
c      write(6,*) "Enables evolution operator computation
      write(6,*) "- EnableEvolutionOperator(bool)"
c      write(6,*) "Lock internal subgrids
      write(6,*) "- LockGrids(bool)"
c      write(6,*) "Switch to the time-like evolution
      write(6,*) "- SetTimeLikeEvolution(bool)"
c      write(6,*) "Switch to the fast evolution
      write(6,*) "- SetFastEvolution(bool)"
c      write(6,*) "Enables the running of the MSbar masses
      write(6,*) "- EnableMassRunning(bool)"
c      write(6,*) "Switch on the small-x resummation
      write(6,*) "- SetSmallxResummation(bool, string la)"
c      write(6,*) "Sets the reference values of $alpha_{s}$ at the reference scale
      write(6,*) "- SetAlphaQCDRef(double alpharef, double Qref)"
c      write(6,*) "Sets the reference values of $alpha$ at the reference scale
      write(6,*) "- SetAlphaQEDRef(double alpharef, double Qref)"
c      write(6,*) "Sets the solution of the beta function
      write(6,*) "- SetAlphaEvolution(string evol)"
c      write(6,*) "Sets the value of LambdaQCD for "nref" flavours
      write(6,*) "- SetLambdaQCDRef(double lambdaref, int nref)"
c      write(6,*) "Sets the solution of the DGLAP equation
      write(6,*) "- SetPDFEvolution(string evolp)"
c      write(6,*) "Sets the minimimum and the maximum energy allowed for the evolution
      write(6,*) "- SetQLimits(double Qmin, double Qmax)"
c      write(6,*) "Sets the FFNS as a default
      write(6,*) "- SetFFNS(int nfl)"
c      write(6,*) "Sets the parameters of the i-th x-space grid
      write(6,*) "- SetGridParameters(int i, int np, int deg, double x)"
c      write(6,*) "Sets the user given i-th x-space grid
      write(6,*) "- SetExternalGrid(int i, int np, int deg, double *x)"
c      write(6,*) "Sets the maximum number of flavours that the evolution 
c      write(6,*) "of alphaQCD and alphaQED can reach
      write(6,*) "- SetMaxFlavourAlpha(int nf)"
c      write(6,*) "Sets the maximum number of flavours that the evolution 
c      write(6,*) "of PDFs can reach.
      write(6,*) "- SetMaxFlavourPDFs(int nf)"
c      write(6,*) "Sets as a default the heavy quark MSbar masses
      write(6,*) "- SetMSbarMasses(double mc, double mb, double mt)"
c      write(6,*) "Sets the number of x-space grids that will be used in the computation
      write(6,*) "- SetNumberOfGrids(int n)"
c      write(6,*) "Sets the name of the PDF set to be used at the initial scale
      write(6,*) "- SetPDFSet(string name)"
c      write(6,*) "Sets the perturbative order of the evolution
      write(6,*) "- SetPerturbativeOrder(int pto)"
c      write(6,*) "Sets as a default the heavy quark pole masses
      write(6,*) "- SetPoleMasses(double mc, double mb, double mt)"
c      write(6,*) "Sets the ratio between renormalization and factorization scales.
      write(6,*) "- SetRenFacRatio(double ratio)"
c      write(6,*) "Sets the replica to be used as initial PDFs (only with a LHAPDF grid)
      write(6,*) "- SetReplica(int nr)"
c      write(6,*) "Sets the FFNS as a default
      write(6,*) "- SetTheory(string theory)"
c      write(6,*) "Sets the VFNS as a default
      write(6,*) "- SetVFNS()"
c      write(6,*) "- /*
c      write(6,*) "- * DIS module
c      write(6,*) "- */
c      write(6,*) "DIS observables
      write(6,*) "- DIS_xsec(double x, double qi, double qf, double y,",
     1           " double pol,"
      write(6,*) "           string proc,string scheme,"
      write(6,*) "           int pto, string pdfset, int irep,"
      write(6,*) "           string target, string proj,"
      write(6,*) "           double *F2, double *F3, double *FL,",
     1           " double *sigma)"
c      write(6,*) "- /*
c      write(6,*) "- * New DIS module
c      write(6,*) "- */
c      write(6,*) "Initialize the new DIS module
      write(6,*) "- InitializeAPFEL_DIS()"
c      write(6,*) "Precompute the structure functions
      write(6,*) "- ComputeStructureFunctionsAPFEL(double Q0, double Q)"
c      write(6,*) "Set the mass scheme for the structure functions
      write(6,*) "- SetMassScheme(string ms)"
c      write(6,*) "Set the polarization
      write(6,*) "- SetPolarizationDIS(double pol)"
c      write(6,*) "- // Set the process of the structure functions (EM, NC or CC)
      write(6,*) "- SetProcessDIS(string pr)"
c      write(6,*) "Set the projectile
      write(6,*) "- SetProjectileDIS(string lept)"
c      write(6,*) "Set the target
      write(6,*) "- SetTargetDIS(string tar)"
c      write(6,*) "Returns the DIS operator times the evolution factors on the grid
      write(6,*) "- ExternalDISOperator(string SF, int ihq,",
     1           " int i, double x, int beta)"
c      write(6,*) "Structure functions
      write(6,*) "- F2light(double x)"
      write(6,*) "- F2charm(double x)"
      write(6,*) "- F2bottom(double x)"
      write(6,*) "- F2top(double x)"
      write(6,*) "- F2total(double x)"
      write(6,*) "- FLlight(double x)"
      write(6,*) "- FLcharm(double x)"
      write(6,*) "- FLbottom(double x)"
      write(6,*) "- FLtop(double x)"
      write(6,*) "- FLtotal(double x)"
      write(6,*) "- F3light(double x)"
      write(6,*) "- F3charm(double x)"
      write(6,*) "- F3bottom(double x)"
      write(6,*) "- F3top(double x)"
      write(6,*) "- F3total(double x)"
c      write(6,*) "Set the value of the Z mass in GeV
      write(6,*) "- SetZMass(double massz)"
c      write(6,*) "Set the value of the W mass in GeV
      write(6,*) "- SetWMass(double massw)"
c      write(6,*) "Set the value of the proton mass in GeV
      write(6,*) "- SetProtonMass(double massp)"
c      write(6,*) "Set the value of sin(theta_W)
      write(6,*) "- SetSinThetaW(double sw)"
c      write(6,*) "Set the absolute value of the entries of the CKM matrix
      write(6,*) "- SetCKM(double vud, double vus, double vub,"
      write(6,*) "         double vcd, double vcs, double vcb,"
      write(6,*) "         double vtd, double vts, double vtb)"
c      write(6,*) "Set the value of the Fermi constant
      write(6,*) "- SetGFermi(double gf)"
c      write(6,*) "Emulator of the FKgenerator
      write(6,*) "- FKSimulator(string obs, double x,",
     1           " double q, double y, int i, int beta)"
c      write(6,*) "Set the observable for FKgenerator simulator
      write(6,*) "- SetFKObservable(string obs)"
c      write(6,*) "Observable according to the FKgenerator naming
      write(6,*) "- FKObservables(string obs, double x,",
     1           " double q, double y)"
      write(6,*) "- ListFunctions()"
      write(6,*) "   "
*
      return
      end
