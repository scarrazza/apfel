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
*
      write(6,*) achar(27)//"[31m---- Functions of the evolution",
     1           " module ----"
      write(6,*) achar(27)//"[0m"
*
      write(6,*) achar(27)//"[33mInitialization functions:"
      write(6,*) achar(27)//"[0m"
      write(6,*) achar(27)//"[34m- InitializeAPFEL():"//achar(27)//"[0m"
      write(6,*) "    initializes the APFEL library. If no settings has"
      write(6,*) "    been specified, it uses the default ones."
      write(6,*) achar(27)//"[34m- EvolveAPFEL(double Q0, double Q):"//
     1           achar(27)//"[0m"
      write(6,*) "    evolves PDFs on the grid to the scale 'Q' [GeV]"
      write(6,*) "    starting from the scale 'Q0' [GeV]."
      write(6,*) achar(27)//"[34m- DeriveAPFEL(double Q):"//
     1           achar(27)//"[0m"
      write(6,*) "    computes the logarithmic derivative with respect"
      write(6,*) "    of 'Q' of PDFs at the scale 'Q' [GeV]."

      write(6,*) achar(27)//"[34m- CachePDFsAPFEL(double Q0):"//
     1           achar(27)//"[0m"
      write(6,*) "    evolves PDFs and cache them on an (x,Q2)-grid"
      write(6,*) "    starting from the scale 'Q0' [GeV]."
      write(6,*) "  "
*
      write(6,*) achar(27)//"[33mSetting functions:"
      write(6,*) achar(27)//"[0m"
      write(6,*) achar(27)//"[34m- SetPerturbativeOrder(int pto):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the perturbative order of the evolution"
      write(6,*) "    (pto = 0,1,2, default 2)."
      write(6,*) achar(27)//"[34m- SetTheory(string theory):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the theory to be used in the evolution"
      write(6,*) "    (theory = 'QCD', 'QUniD', default 'QCD')"
      write(6,*) achar(27)//"[34m- EnableNLOQEDCorrections(bool):"//
     1           achar(27)//"[0m"
      write(6,*) "    enables/disables the NLO QED corrections. If"
      write(6,*) "    disabled only LO QED are included. This is active"
      write(6,*) "    only for the 'QUniD' solution (default false)."
      write(6,*) achar(27)//"[34m- SetVFNS():"//achar(27)//"[0m"
      write(6,*) "    sets the Variable-Flavour Number Scheme."
      write(6,*) achar(27)//"[34m- SetFFNS(int nfl):"//achar(27)//"[0m"
      write(6,*) "    sets the Fixed-Flavour Number Scheme with 'nfl'"
      write(6,*) "    active flavours."
      write(6,*) achar(27)//"[34m- SetAlphaQCDRef(double alpharef,",
     1           " double Qref):"//achar(27)//"[0m"
      write(6,*) "    sets the reference values of alphas at the scale"
      write(6,*) "    'Qref' [GeV] to 'alpharef' (default 'alpharef' ="
      write(6,*) "    0.35, 'Qref' = sqrt(2) GeV)."
      write(6,*) achar(27)//"[34m- SetAlphaQEDRef(double alpharef,",
     1           " double Qref):"//achar(27)//"[0m"
      write(6,*) "    sets the reference values of alpha at the scale"
      write(6,*) "    'Qref' [GeV] to 'alpharef' (default 'alpharef' ="
      write(6,*) "    7.496252d-3, 'Qref' = 1.777 GeV)."
      write(6,*) achar(27)//"[34m- SetLambdaQCDRef(double lambdaref,",
     1           " int nref):"//achar(27)//"[0m"
      write(6,*) "    sets the value of LambdaQCD [GeV] with 'nref'"
      write(6,*) "    flavours to 'lambdaref' (default 'lambdaref'"
      write(6,*) "    = 0.220, 'nref' = 5)"
      write(6,*) achar(27)//"[34m- SetPoleMasses(double mc, double mb,",
     1           " double mt):"//achar(27)//"[0m"
      write(6,*) "    sets the values of the heavy quark thresholds"
      write(6,*) "    in GeV in the Pole-mass scheme (default 'mc' "
      write(6,*) "    = sqrt(2) GeV, 'mb' = 4.5 GeV, 'mt' = 175 GeV)."
      write(6,*) achar(27)//"[34m- SetMSbarMasses(double mc,",
     1           " double mb, double mt):"//achar(27)//"[0m"
      write(6,*) "    sets the values of the heavy quark thresholds"
      write(6,*) "    in GeV in the MSbar scheme."
      write(6,*) achar(27)//"[34m- SetMassScaleReference(double Qc,",
     1           " double Qb, double Qt):"//achar(27)//"[0m"
      write(6,*) "    sets the reference scales in GeV at which heavy"
      write(6,*) "    quark masses are given. This has no effect if the"
      write(6,*) "    pole masses are used."
      write(6,*) achar(27)//"[34m- SetMassMatchingScales(double kmc,",
     1           " double kmb, double kmt):"//achar(27)//"[0m"
      write(6,*) "    sets the ratios between heavy quark masses and"
      write(6,*) "    heavy quark matching thresholds"
      write(6,*) "    (default 'kmc' = 1, 'kmb' = 1, 'kmt' = 1)."
      write(6,*) achar(27)//"[34m- SetTauMass(double mtau):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the values of the tau lepton in GeV"
      write(6,*) "    (default 'mtau' = 1.777 GeV)"
      write(6,*) achar(27)//"[34m- EnableMassRunning(bool):"//
     1           achar(27)//"[0m"
      write(6,*) "    enables/disables the running of the MSbar masses"
      write(6,*) "    (default true)."
      write(6,*) achar(27)//"[34m- SetMaxFlavourAlpha(int nf):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the maximum number of active flavours in the"
      write(6,*) "    couplings evolution (including the masses)"
      write(6,*) "    (default 'nf' = 6)."
      write(6,*) achar(27)//"[34m- SetMaxFlavourPDFs(int nf):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the maximum number of active flavours in the"
      write(6,*) "    PDF evolution (default 'nf' = 6)."
      write(6,*) achar(27)//"[34m- SetRenFacRatio(double ratio):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the ratio between factorization and"
      write(6,*) "    renormalization scales both in GeV (default"
      write(6,*) "    'ratio' = 1)."
      write(6,*) achar(27)//"[34m- SetTimeLikeEvolution(bool):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the time-like evolution (frag. functions)"
      write(6,*) "    (default false)."
      write(6,*) achar(27)//"[34m- SetPolarizedEvolution(bool):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the polarized evolution (default false)"
      write(6,*) achar(27)//"[34m- SetSmallxResummation(bool,",
     1           " string la):"//achar(27)//"[0m"
      write(6,*) "    sets the the small-x resummation at 'la' log"
      write(6,*) "    accuracy ('la' = 'LL', 'NLL') (default false)."
      write(6,*) achar(27)//"[34m- SetAlphaEvolution(string evol):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the solution of the beta function equations"
      write(6,*) "    for the running couplings ('evol' = 'exact',"
      write(6,*) "    'expanded', 'lambda') (default 'evol' = 'exact')"
      write(6,*) achar(27)//"[34m- SetPDFEvolution(string evolp):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the solution of the DGLAP equations for PDFs"
      write(6,*) "    ('evolp' = 'exactmu', 'exactalpha', 'expandalpha'"
      write(6,*) "    'truncated') (default 'evolp' = 'exactmu')"
      write(6,*) achar(27)//"[34m- SetEpsilonTruncation(double eps):"//
     1           achar(27)//"[0m"
      write(6,*) "    if the truncated evolution for PDFs has been"
      write(6,*) "    chosen, it sets the truncation parameter 'eps'."
      write(6,*) achar(27)//"[34m- SetPDFSet(string name):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the PDF set to be evolved. 'name' can be an"
      write(6,*) "    LHAPDF set (must finish with '.LHgrid'), some of"
      write(6,*) "    the internal set or an external set (default"
      write(6,*) "    'name' = 'ToyLH')."
      write(6,*) achar(27)//"[34m- SetReplica(int nr):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the replica/member of the LHAPDF set to be"
      write(6,*) "    evolved (default 'nr' = 0)."
      write(6,*) achar(27)//"[34m- SetQLimits(double Qmin,",
     1           " double Qmax):"//achar(27)//"[0m"
      write(6,*) "    sets the range where it is possible to perform"
      write(6,*) "    the evolution (default 'Qmin' = 0.5 GeV, 'Qmax'"
      write(6,*) "    = 1000 GeV)."
      write(6,*) achar(27)//"[34m- SetNumberOfGrids(int n):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the number of subgrids 'n' (default 3)"
      write(6,*) achar(27)//"[34m- SetGridParameters(int i, int n,",
     1           " int deg, double x):"//achar(27)//"[0m"
      write(6,*) "    sets the parameters of the i-th subgrid. 'n' ="
      write(6,*) "    number of intevals, 'deg' = interpolation order,"
      write(6,*) "    'x' lower bound of the grid (the upper bound is"
      write(6,*) "    always 1)."
      write(6,*) achar(27)//"[34m- SetQGridParameters(int nQ,",
     1           " int degQ):"//achar(27)//"[0m"
      write(6,*) "    sets the parameters of the Q-grid. 'nQ' = number"
      write(6,*) "    of intevals and 'degQ' = interpolation order."
      write(6,*) "    (default: 'nQ' = 100, 'degQ' = 3, relevant only"
      write(6,*) "    for the cached evolution)."
      write(6,*) achar(27)//"[34m- SetExternalGrid(int i, int np,",
     1           " int deg, double *x):"//achar(27)//"[0m"
      write(6,*) "    sets the external grid in the position 'i' with"
      write(6,*) "    'np' intervals, interpolation degree 'deg'. 'x'"
      write(6,*) "    must be a one-dimentional array with upper bound"
      write(6,*) "    in 1 (there cannot be more than 1 external grid)."
      write(6,*) achar(27)//"[34m- SetFastEvolution(bool):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the fast PDF evolution (default true)."
      write(6,*) achar(27)//"[34m- GetVersion():"//achar(27)//"[0m"
      write(6,*) "    returns the APFEL version in use."
      write(6,*) achar(27)//"[34m- EnableWelcomeMessage(bool):"//
     1           achar(27)//"[0m"
      write(6,*) "    enables the printing of the welcome message with"
      write(6,*) "    the APFEL banner and the report of the evolution"
      write(6,*) "    parameters (default true)."
      write(6,*) achar(27)//"[34m- EnableEvolutionOperator(bool):"//
     1           achar(27)//"[0m"
      write(6,*) "    enables the computation of the external evolution"
      write(6,*) "    parameters (default false)."
      write(6,*) achar(27)//"[34m- EnableLeptonEvolution(bool):"//
     1           achar(27)//"[0m"
      write(6,*) "    enables the evolution of the lepton PDFs when the"
      write(6,*) "    fast QUniD is used (default false)."
      write(6,*) achar(27)//"[34m- LockGrids(bool):"//achar(27)//"[0m"
      write(6,*) "    locks the subgrids (default false)."
      write(6,*) achar(27)//"[34m- CleanUp():"//achar(27)//"[0m"
      write(6,*) "    resets all the evolution parameters to the"
      write(6,*) "    default settings."
      write(6,*) achar(27)//"[34m- SetLHgridParameters(int nx,",
     1           " int nxm, double xmin, double xm, double xmax,",
     2           " int nq2,"," double q2min, double q2max):"//
     3           achar(27)//"[0m"
      write(6,*) "    sets the parameters of the grid over which PDFs"
      write(6,*) "    will be tabulated in the LHAPDF format."
      write(6,*) achar(27)//"[34m- ListFunctions():"//achar(27)//"[0m"
      write(6,*) "    lists all the functions available in APFEL."
      write(6,*) "   "
*
      write(6,*) achar(27)//"[33mOutput functions:"//achar(27)//"[0m"
      write(6,*) achar(27)//"[0m"
      write(6,*) achar(27)//"[34m- xPDF(int i, double x) and ",
     1           "xgamma(double x):"//achar(27)//"[0m"
      write(6,*) "    return 'x' times the i-th and the photon PDF"
      write(6,*) "    in 'x' at the final scale 'Q' [GeV] defined in"
      write(6,*) "    'EvolveAPFEL'."
      write(6,*) achar(27)//"[34m- xPDFall(double x, double *xf):"//
     1           achar(27)//"[0m"
      write(6,*) "    returns at once 'x' times all the PDF in the"
      write(6,*) "    array xf[-6:6] computed in 'x' at the final scale"
      write(6,*) "    'Q' [GeV] defined in 'EvolveAPFEL'."
      write(6,*) achar(27)//"[34m- xPDFallPhoton(double x,",
     1           " double *xf):"//achar(27)//"[0m"
      write(6,*) "    returns at once 'x' times all the PDF, including"
      write(6,*) "    the photon, in the array xf[-6:7] computed in 'x'"
      write(6,*) "    at the final scale 'Q' [GeV] defined in "
      write(6,*) "    'EvolveAPFEL'."
      write(6,*) achar(27)//"[34m- xPDFj(int i, double x) and",
     1           " xgammaj(double x):"//achar(27)//"[0m"
      write(6,*) "    return 'x' times the i-th and the photon PDF"
      write(6,*) "    in 'x' at the final scale 'Q' [GeV] defined in"
      write(6,*) "    'EvolveAPFEL' interpolated on the joint grid."
      write(6,*) achar(27)//"[34m- dxPDF(int i, double x) and",
     1           " dxgamma(double x):"//achar(27)//"[0m"
      write(6,*) "    return 'x' times the derivative in ln(Q^2) of"
      write(6,*) "    the i-th and the photon PDF in 'x' at the scale"
      write(6,*) "    'Q' [GeV] defined in 'DeriveAPFEL'."
      write(6,*) achar(27)//"[34m- NPDF(int i, int N)",
     1           " and Ngamma(int N):"//achar(27)//"[0m"
      write(6,*) "    return the N-th moment of the i-th and the"
      write(6,*) "    photon PDF the final scale 'Q' [GeV] defined in"
      write(6,*) "    'EvolveAPFEL'."
      write(6,*) achar(27)//"[34m- xLepton(int i, double x):"//
     1           achar(27)//"[0m"
      write(6,*) "    return 'x' times the i-th lepton PDF in 'x'"
      write(6,*) "    at the final scale 'Q' [GeV] defined in"
      write(6,*) "    'EvolveAPFEL'."
      write(6,*) achar(27)//"[34m- xLeptonj(int i, double x):"//
     1           achar(27)//"[0m"
      write(6,*) "    return 'x' times the i-th lepton PDF in 'x'"
      write(6,*) "    at the final scale 'Q' [GeV] defined in"
      write(6,*) "    'EvolveAPFEL' interpolated on the joint grid."
      write(6,*) achar(27)//"[34m- xPDFxQ(int i, double x, double Q):"//
     1           achar(27)//"[0m"
      write(6,*) "    return 'x' times the i-th PDF (inclunding quarks,"
      write(6,*) "    gluon, photon and leptons) in 'x' at the final"
      write(6,*) "    scale 'Q' [GeV]. This function requires that"
      write(6,*) "    'CachePDFsAPFEL' to be called in advance."
      write(6,*) achar(27)//"[34m- xPDFxQall(double x, double Q,",
     1           " double *xf):"//achar(27)//"[0m"
      write(6,*) "    returns at once 'x' times all the PDF in the"
      write(6,*) "    array xf[-6:6] computed in 'x' at the scale"
      write(6,*) "    'Q' [GeV]. This function requires that"
      write(6,*) "    'CachePDFsAPFEL' to be called in advance."
      write(6,*) achar(27)//"[34m- LUMI(int i, int j, double S):"//
     1           achar(27)//"[0m"
      write(6,*) "    returns the partonic luminosity of the i-th and"
      write(6,*) "    j-th partons for the CoM energy S [GeV^2] for the"
      write(6,*) "    final invariant mass Mx = Q [GeV] defined in "
      write(6,*) "    'EvolveAPFEL'."
      write(6,*) achar(27)//"[34m- AlphaQCD(double Q):"//
     1           achar(27)//"[0m"
      write(6,*) "    returns the QCD strong coupling alpha_s at the"
      write(6,*) "    scale 'Q' [GeV]."
      write(6,*) achar(27)//"[34m- AlphaQED(double Q):"//
     1           achar(27)//"[0m"
      write(6,*) "    returns the QED coupling alpha at the scale"
      write(6,*) "    'Q' [GeV]."
      write(6,*) achar(27)//"[34m- HeavyQuarkMass(int i,double Q):"//
     1           achar(27)//"[0m"
      write(6,*) "    returns the mass of the i-th heavy quark"
      write(6,*) "    (i = 4,5,6) scale 'Q' [GeV] (the masses run only"
      write(6,*) "    when using the MSbar scheme)."
      write(6,*) achar(27)//"[34m- GetLambdaQCD(int i):"//
     1           achar(27)//"[0m"
      write(6,*) "    returns the value of LambdaQCD with in GeV with"
      write(6,*) "    'i' active flavours (i = 4,5,6)."
      write(6,*) achar(27)//"[34m- GetThreshold(int i):"//
     1           achar(27)//"[0m"
      write(6,*) "    returns the value of the i-th heavy quark"
      write(6,*) "    threhold in GeV (i = 4,5,6)."
      write(6,*) achar(27)//"[34m- GetMaxFlavourAlpha():"//
     1           achar(27)//"[0m"
      write(6,*) "    return the maximum number of active flavours in"
      write(6,*) "    the couplings evolution."
      write(6,*) achar(27)//"[34m- GetMaxFlavourPDFs():"//
     1           achar(27)//"[0m"
      write(6,*) "    return the maximum number of active flavours in"
      write(6,*) "    the PDF evolution."
      write(6,*) achar(27)//"[34m- nIntervals():"//achar(27)//"[0m"
      write(6,*) "    returns the number of intervals of the joint"
      write(6,*) "    grid."
      write(6,*) achar(27)//"[34m- xGrid(int alpha):"//achar(27)//"[0m"
      write(6,*) "    returns the value of 'x' on the alpha-th node of"
      write(6,*) "    the joint grid."
      write(6,*) achar(27)//"[34m- GetPerturbativeOrder():"//
     1           achar(27)//"[0m"
      write(6,*) "    returns the perturbative order set for the"
      write(6,*) "    evolution."
      write(6,*) achar(27)//"[34m- GetMuF0():"//
     1           achar(27)//"[0m"
      write(6,*) "    returns the initial factorizationn scale."
      write(6,*) achar(27)//"[34m- GetMuF():"//
     1           achar(27)//"[0m"
      write(6,*) "    returns the final factorizationn scale."
      write(6,*) achar(27)//"[34m- ExternalEvolutionOperator(",
     1           "string fname, int i, int j, double x, int beta):"//
     2           achar(27)//"[0m"
      write(6,*) "    returns the PDF evolution operator."
      write(6,*) achar(27)//"[34m- ExternalEvolutionMatrixEv2Ev(",
     1           "int i, int j, int alpha, int beta):"//
     2           achar(27)//"[0m"
      write(6,*) "    returns the PDF evolution matrix from"
      write(6,*) "    evolution to evolution basis."
      write(6,*) achar(27)//"[34m- ExternalEvolutionMatrixEv2Ph(",
     1           "int i, int j, int alpha, int beta):"//
     2           achar(27)//"[0m"
      write(6,*) "    returns the PDF evolution matrix from"
      write(6,*) "    evolution to physical basis."
      write(6,*) achar(27)//"[34m- ExternalEvolutionMatrixPh2Ph(",
     1           "int i, int j, int alpha, int beta):"//
     2           achar(27)//"[0m"
      write(6,*) "    returns the PDF evolution matrix from"
      write(6,*) "    physical to physical basis."
      write(6,*) achar(27)//"[34m- ExternalSplittingFunctions(",
     1           "string fname, int pt, int nf, int i, int j,",
     2           " double x, int beta):"//achar(27)//"[0m"
      write(6,*) "    returns the QCD splitting functions on the"
      write(6,*) "    joint interpolation grid."
      write(6,*) achar(27)//"[34m- LHAPDFgrid(int Nrep, double Qin,",
     1           " string fname):"//achar(27)//"[0m"
      write(6,*) "    produces a PDF interpolation grid in the LHAPDF"
      write(6,*) "    format."
      write(6,*) achar(27)//"[34m- LHAPDFgridDerivative(int Nrep,",
     1           " string fname):"//achar(27)//"[0m"
      write(6,*) "    produces an interpolation grid in the LHAPDF"
      write(6,*) "    format for the derived PDFs."
      write(6,*) "  "
*
      write(6,*) achar(27)//"[31m---- Functions of the DIS module ----"
      write(6,*) achar(27)//"[0m"
*
      write(6,*) achar(27)//"[33mInitialization functions:"
      write(6,*) achar(27)//"[0m"
      write(6,*) achar(27)//"[34m- InitializeAPFEL_DIS():"//
     1           achar(27)//"[0m"
      write(6,*) "    initializes the DIS module. If no settings has"
      write(6,*) "    been specified, it uses the default ones."
      write(6,*) achar(27)//"[34m- ComputeStructureFunctionsAPFEL(",
     1           "double Q0, double Q):"//achar(27)//"[0m"
      write(6,*) "    computes the DIS structure functions on the grid"
      write(6,*) "    at the scale 'Q' [GeV] applying also the PDF"
      write(6,*) "    evolution from the initial scale 'Q0' [GeV]."
      write(6,*) "  "
      write(6,*) achar(27)//"[33mSetting functions:"
      write(6,*) achar(27)//"[0m"
      write(6,*) achar(27)//"[34m- SetMassScheme(string ms):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the mass scheme to be used to compute the"
      write(6,*) "    structure functions ('ms' = 'ZM-VFNS', 'FFNS',"
      write(6,*) "    'FONLL-A', 'FONLL-B', 'FONLL-C', default 'ms' ="
      write(6,*) "    'ZM-VFNS')."
      write(6,*) achar(27)//"[34m- SetPolarizationDIS(double pol):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the beam polarization (default 'pol' = 0)."
      write(6,*) achar(27)//"[34m- SetProcessDIS(string pr):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets process ('pr' = 'EM', 'NC', 'CC', default"
      write(6,*) "    'pr' = 'EM')."
      write(6,*) achar(27)//"[34m- SetNCComponent(string cm):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the component of the NC structure functions"
      write(6,*) "    ('cm' = 'gg', 'gZ', 'ZZ', 'al', default"
      write(6,*) "    'cm' = 'al')."
      write(6,*) achar(27)//"[34m- SetProjectileDIS(string lept):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the projectile ('lept' = 'electron',"
      write(6,*) "    'positron', 'neutrino', 'antineutrino', default"
      write(6,*) "    'lept' = 'electron')."
      write(6,*) achar(27)//"[34m- SetTargetDIS(string tar):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the target ('tar' = 'proton', 'neutron',"
      write(6,*) "    'isoscalar', 'iron', default 'tar' = 'proton')"
      write(6,*) achar(27)//"[34m- SetZMass(double massz):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the value of the mass of the Z boson"
      write(6,*) "    (default 'massz' = 91.1876 GeV)."
      write(6,*) achar(27)//"[34m- SetWMass(double massw):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the value of the mass of the W boson"
      write(6,*) "    (default 'massw' = 80.385 GeV)."
      write(6,*) achar(27)//"[34m- SetProtonMass(double massp):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the value of the mass of the proton"
      write(6,*) "    (default 'massp' = 0.938272046 GeV)."
      write(6,*) achar(27)//"[34m- SetSin2ThetaW(double sw):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the value of sin^2(theta_W)"
      write(6,*) "    (default 'sw' = 0.23126)."
      write(6,*) achar(27)//"[34m- SetGFermi(double gf):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the value of Fermi constant"
      write(6,*) "    (default 'gf' = 1.1663787e-5)."
      write(6,*) achar(27)//"[34m- SetCKM(double vud, double vus,",
     1           " double vub,"
      write(6,*) "         double vcd, double vcs, double vcb,"
      write(6,*) "         double vtd, double vts, double vtb):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the absolute value of the entries of the"
      write(6,*) "    CKM matrix"
      write(6,*) "    (default: 0.97427d0, 0.22536d0, 0.00355d0,"
      write(6,*) "              0.22522d0, 0.97343d0, 0.04140d0,"
      write(6,*) "              0.00886d0, 0.04050d0, 0.99914d0)."
      write(6,*) achar(27)//"[34m- SetRenQRatio(double ratio):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the ratio muR / Q (default 1)"
      write(6,*) achar(27)//"[34m- SetFacQRatio(double ratio):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the ratio muF / Q (default 1)"
      write(6,*) achar(27)//"[34m- SetScaleVariationProcedure(int",
     1           " scp):"//achar(27)//"[0m"
      write(6,*) "    sets the scale variation procedure (default: 0)."
      write(6,*) "    The options are:"
      write(6,*) "    - 0: scale variations in evolution"
      write(6,*) "         and structure functions."
      write(6,*) "    - 1: scale variations in structure functions"
      write(6,*) "         only."
      write(6,*) achar(27)//"[34m- EnableDynamicalScaleVariations(",
     1           "bool):"//achar(27)//"[0m"
      write(6,*) "    enables or disables the possibility to perform"
      write(6,*) "    fact/ren scale variations point by point without"
      write(6,*) "    requiring the ratio \mu_{R,F} / Q to be constant."
      write(6,*) "    Limitations: \mu_F = \mu_R and slower code."
      write(6,*) achar(27)//"[34m- EnableIntrinsicCharm(",
     1           "bool):"//achar(27)//"[0m"
      write(6,*) "    enables or disables the intrinsic charm"
      write(6,*) "    contributions to the massive structure functions."      
      write(6,*) achar(27)//"[34m- EnableTargetMassCorrections(bool):"//
     1           achar(27)//"[0m"
      write(6,*) "    enables or disables the target mass corrections"
      write(6,*) "    to the DIS structure functions due to the finite"
      write(6,*) "    mass of the proton."
      write(6,*) achar(27)//"[34m- EnableDampingFONLL(bool):"//
     1           achar(27)//"[0m"
      write(6,*) "    enables or disables the damping factor when the"
      write(6,*) "    FONLL structure functions are computed."
      write(6,*) achar(27)//"[34m- SetDampingPowerFONLL(int dpc,",
     1     " int dpb, int dpt):"//achar(27)//"[0m"
      write(6,*) "    set the power with which damping factor"
      write(6,*) "    suppresses the subleading terms for charm,"
      write(6,*) "    bottom, and top, separately."
      write(6,*) achar(27)//"[34m- SelectCharge(string selch):"//
     1           achar(27)//"[0m"
      write(6,*) "    selects one particular charge in the NC structure"
      write(6,*) "    functions ('selch' = 'down', 'up', 'strange',"
      write(6,*) "    'charm', 'bottom', 'top', 'all', default "
      write(6,*) "    'selch' = 'all')"
      write(6,*) achar(27)//"[34m- SetPropagatorCorrection(double dr):"
     1           //achar(27)//"[0m"
      write(6,*) "    sets the correction to the Z propagator involved"
      write(6,*) "    in the NC DIS structure functions"
      write(6,*) "    (default 'dr' = 0)."
      write(6,*) achar(27)//"[34m- SetEWCouplings(double vd,",
     1           " double vu, double ad, double au):"//
     1           achar(27)//"[0m"
      write(6,*) "    sets the vector and axial couplings of the up-"
      write(6,*) "    and down-type quarks. If they are not set by the"
      write(6,*) "    user the standard couplinglings are used."
      write(6,*) "  "
*
      write(6,*) achar(27)//"[33mOutput functions:"
      write(6,*) achar(27)//"[0m"
      write(6,*) achar(27)//"[34m- F2light(double x),",
     1           " F2charm(double x), F2bottom(double x),"
      write(6,*) "  F2top(double x), F2total(double x):"//
     1           achar(27)//"[0m"
      write(6,*) achar(27)//"[34m- FLlight(double x),",
     1           " FLcharm(double x), FLbottom(double x),"
      write(6,*) "  FLtop(double x), FLtotal(double x):"//
     1           achar(27)//"[0m"
      write(6,*) achar(27)//"[34m- F3light(double x),",
     1           " F3charm(double x), F3bottom(double x),"
      write(6,*) "  F3top(double x), F3total(double x):"//
     1           achar(27)//"[0m"
      write(6,*) "    return the F2, FL and xF3 structure functions in"
      write(6,*) "    'x' at the final scale 'Q' [GeV] defined in"
      write(6,*) "    'ComputeStructureFunctionsAPFEL'."
      write(6,*) achar(27)//"[34m- g1light(double x),",
     1           " g1charm(double x), g1bottom(double x),"
      write(6,*) "  g1top(double x), g1total(double x):"//
     1           achar(27)//"[0m"
      write(6,*) achar(27)//"[34m- gLlight(double x),",
     1           " gLcharm(double x), gLbottom(double x),"
      write(6,*) "  gLtop(double x), gLtotal(double x):"//
     1           achar(27)//"[0m"
      write(6,*) achar(27)//"[34m- g4light(double x),",
     1           " g4charm(double x), g4bottom(double x),"
      write(6,*) "  g4top(double x), g4total(double x):"//
     1           achar(27)//"[0m"
      write(6,*) "    return the g1, gL and xg4 polarised structure"
      write(6,*) "    functions in 'x' at the final scale 'Q' [GeV]"
      write(6,*) "     defined in 'ComputeStructureFunctionsAPFEL'."
      write(6,*) achar(27)//"[34m- GetZMass():"//achar(27)//"[0m"
      write(6,*) "    returns the value of the mass of the Z boson"
      write(6,*) achar(27)//"[34m- GetWMass():"//achar(27)//"[0m"
      write(6,*) "    returns the value of the mass of the W boson"
      write(6,*) achar(27)//"[34m- GetProtonMass():"//achar(27)//"[0m"
      write(6,*) "    returns the value of the mass of the proton"
      write(6,*) achar(27)//"[34m- GetSin2ThetaW():"//achar(27)//"[0m"
      write(6,*) "    returns the value of sin^2(theta_W)"
      write(6,*) achar(27)//"[34m- GetGFermi():"//achar(27)//"[0m"
      write(6,*) "    returns the value of Fermi constant"
      write(6,*) achar(27)//"[34m- GetCKM(int u, int d):"//
     1           achar(27)//"[0m"
      write(6,*) "    returns the absolute value of the (u,d) entry"
      write(6,*) "    of the CKM matrix"
      write(6,*) achar(27)//"[34m- GetSIATotalCrossSection(int pto,",
     1           " double Q):"//achar(27)//"[0m"
      write(6,*) "    returns the SIA total cross section in natural"
      write(6,*) "    units at the perturbative order 'pto' and at the"
      write(6,*) "    scale 'Q' in GeV (only for time-like evolution)."
      write(6,*) achar(27)//"[34m- ExternalDISOperator(string SF,",
     1           " int ihq, int i, double x, int beta):"//
     2           achar(27)//"[0m"
      write(6,*) "    returns the DIS operators."
      write(6,*) "  "
*
      write(6,*) achar(27)//
     1           "[33mSpecial functions for the production of FK",
     2           " tables:"//achar(27)//"[0m"
      write(6,*) "  "
      write(6,*) achar(27)//"[34m- SetFKObservable(string obs)"
      write(6,*) achar(27)//"[34m- GetFKObservable()"
      write(6,*) achar(27)//"[34m- FKSimulator(string obs, double x,",
     1           " double q, double y, int i, int beta)"
      write(6,*) achar(27)//"[34m- FKObservables(double x,",
     1           " double q, double y)"
      write(6,*) achar(27)//"[34m- ComputeHardCrossSectionsDY(",
     1           "string inputfile, string outputfile)"
      write(6,*) achar(27)//"[34m- ComputeFKTables(string inputfile, ",
     1           "double Q0, int flmap[196])"//achar(27)//"[0m"
      
      write(6,*) "   "
*
      return
      end
