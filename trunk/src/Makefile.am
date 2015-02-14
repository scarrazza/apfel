noinst_LTLIBRARIES = libAPFELCore.la

INIT = 	Evolution/apfelversion.cc Evolution/WelcomeMessage.f Evolution/initParameters.f \
	Evolution/initGrid.f Evolution/initGridAlpha.f Evolution/initIntegralsQCD.f \
	Evolution/initIntegralsQCDRes.f Evolution/initIntegralsQED.f Evolution/initIntegralsMatching.f \
	Evolution/initPDFs.f Evolution/InitializeAPFEL.f Evolution/EvolveAPFEL.f Evolution/DeriveAPFEL.f \
	Evolution/initHELL.f

SETTINGS = Evolution/SetAlphaQCDRef.f Evolution/SetAlphaQEDRef.f Evolution/SetQLimits.f \
	Evolution/SetPerturbativeOrder.f Evolution/SetFFNS.f Evolution/SetGridParameters.f Evolution/SetExternalGrid.f \
	Evolution/SetMSbarMasses.f Evolution/SetNumberOfGrids.f Evolution/SetPDFSet.f Evolution/SetPoleMasses.f \
	Evolution/SetRenFacRatio.f Evolution/SetTheory.f Evolution/SetVFNS.f Evolution/SetMaxFlavourPDFs.f \
	Evolution/SetMaxFlavourAlpha.f Evolution/SetReplica.f Evolution/AlphaQCD.f Evolution/AlphaQED.f Evolution/HeavyQuarkMass.f \
	Evolution/LHAPDFgrid.f Evolution/CleanUp.f Evolution/EnableWelcomeMessage.f Evolution/EnableEvolutionOperator.f \
	Evolution/ExternalEvolutionOperator.f Evolution/LockGrids.f Evolution/SetTimeLikeEvolution.f Evolution/SetAlphaEvolution.f \
	Evolution/SetLambdaQCDRef.f Evolution/SetPDFEvolution.f Evolution/hknsff07.f Evolution/fDSS.f Evolution/LHAPDFgridDerivative.f \
	Evolution/SetSmallxResummation.f Evolution/SetFastEvolution.f Evolution/EnableMassRunning.f Evolution/ListFunctions.f

INTER = Evolution/interpolants.f

SPLITTING_FUNCTIONS = Evolution/SplittingFunctions.f Evolution/Hellx.f

COUPLINGS = Evolution/a_QCD.f Evolution/a_QED.f Evolution/MSbarmass.f Evolution/lambda.f

EVOLUTIONQCD = Evolution/EvolutionQCD.f Evolution/EvolutionOperatorsQCD.f Evolution/IdentityOperatorsQCD.f \
	Evolution/EqualOperatorsQCDnf.f Evolution/integrandsQCD.f Evolution/integralsQCD.f Evolution/RSLintegralsQCD.f \
	Evolution/odeintnsQCD.f Evolution/odeintsgQCD.f Evolution/MatchingConditions.f Evolution/integrandsMatching.f \
	Evolution/integralsMatching.f Evolution/RSLintegralsMatching.f Evolution/MatchPDFs.f Evolution/JoinOperatorsQCD.f \
	Evolution/EvolutionUnified.f Evolution/EvolutionOperatorsUnified.f Evolution/IdentityOperatorsUnified.f \
	Evolution/EqualOperatorsUnifiednf.f Evolution/EvolveUni.f Evolution/odeintnsUnified.f Evolution/odeintsgUnified.f \
	Evolution/DerivativeOperatorsQCD.f Evolution/odeintnsQCDf.f Evolution/odeintsgQCDf.f \
	Evolution/odeintnsUnifiedf.f Evolution/odeintsgUnifiedf.f

EVOLUTIONQED = Evolution/EvolutionQED.f Evolution/EvolutionOperatorsQED.f Evolution/IdentityOperatorsQED.f \
	Evolution/EqualOperatorsQEDnf.f Evolution/integrandsQED.f Evolution/integralsQED.f Evolution/RSLintegralsQED.f \
	Evolution/odeintnsQED.f Evolution/odeintsgQED.f Evolution/DerivativeOperatorsQED.f \
	Evolution/odeintnsQEDf.f Evolution/odeintsgQEDf.f

PDF = Evolution/PDFphys2evQCD.f Evolution/PDFevQCD2phys.f Evolution/PDFphys2evQED.f Evolution/PDFevQED2phys.f \
	Evolution/EvolvePDFs.f Evolution/switchGluonPhoton.f Evolution/EvolveQCD.f Evolution/EvolveQED.f Evolution/toyLHPDFs.f \
	Evolution/xPDF.f Evolution/xgamma.f Evolution/NPDF.f Evolution/Ngamma.f Evolution/LUMI.f Evolution/JoinGrids.f \
	Evolution/PDFevUni2evQCD.f Evolution/PDFevQCD2evUni.f Evolution/PDFevUni2phys.f Evolution/PDFphys2evUni.f \
	Evolution/dxPDF.f Evolution/dxgamma.f Evolution/DeriveQCD.f Evolution/DeriveQED.f Evolution/DerivePDFs.f \
	Evolution/xGrid.f Evolution/nIntervals.f Evolution/ExternalSetAPFEL.f

UTILITIES = core/dgauss.f core/wgplg.f core/ddilog.f core/gamma.f core/beta.f

DIS = DIS/cc_dis.f DIS/FF0CoefficientFunctions.f DIS/integrands_cc.f DIS/nc_dis.f DIS/xcoeffunc_cc.f DIS/xpdfs.f DIS/dis_xsec.f DIS/FFCoefficientFunctions.f DIS/integrands_nc.f DIS/PDFparam.f DIS/xcoeffunc_nc.f DIS/ZMCoefficientFunctions.f

DIS_NEW = DIS_new/InitializeAPFEL_DIS.f DIS_new/initParametersDIS.f DIS_new/SetMassScheme.f DIS_new/initIntegralsDIS.f DIS_new/ZeroMassCoefficientFunctions.f DIS_new/RSLintegralsDIS.f DIS_new/ComputeStructureFunctionsAPFEL.f DIS_new/F2total.f DIS_new/F2light.f DIS_new/F2charm.f DIS_new/F2bottom.f DIS_new/F2top.f DIS_new/FLtotal.f DIS_new/FLlight.f DIS_new/FLcharm.f DIS_new/FLbottom.f DIS_new/FLtop.f DIS_new/F3total.f DIS_new/F3light.f DIS_new/F3charm.f DIS_new/F3bottom.f DIS_new/F3top.f DIS_new/integrandsDIS.f DIS_new/SetProcessDIS.f DIS_new/ComputeChargesDIS.f DIS_new/SetPolarizationDIS.f DIS_new/SetProjectileDIS.f DIS_new/SetTargetDIS.f DIS_new/MassiveCoefficientFunctions.f DIS_new/MassiveZeroCoefficientFunctions.f DIS_new/ComputeDISOperators.f DIS_new/ConvolutePDFsWithDISOperators.f DIS_new/JoinDISOperators.f DIS_new/ConvoluteEvolutionWithDISOperators.f DIS_new/ExternalDISOperator.f DIS_new/SetCKM.f DIS_new/SetGFermi.f DIS_new/SetProtonMass.f DIS_new/SetSinThetaW.f DIS_new/SetWMass.f DIS_new/SetZMass.f DIS_new/FKSimulator.f DIS_new/SetFKObservable.f DIS_new/FKObservables.f

EXTRA_DIST = commons

HELL = HELL/src/hell.cc HELL/src/DGLAPevol.cc HELL/src/hell_wrapper.cc

libAPFELCore_la_SOURCES = $(INIT) $(SETTINGS) $(INTER) $(SPLITTING_FUNCTIONS) $(COUPLINGS) $(EVOLUTIONQCD) $(EVOLUTIONQED) $(PDF) $(UTILITIES) $(DIS) $(DIS_NEW) $(HELL)

AM_CXXFLAGS += -I$(shell lhapdf-config --incdir)
