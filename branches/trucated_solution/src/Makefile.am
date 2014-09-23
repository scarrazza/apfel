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
	Evolution/SetSmallxResummation.f Evolution/SetFastEvolution.f Evolution/SetEpsTrunc.f

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
	Evolution/xGrid.f Evolution/nIntervals.f

UTILITIES = core/dgauss.f core/wgplg.f core/ddilog.f

DIS = DIS/cc_dis.f DIS/FF0CoefficientFunctions.f DIS/integrands_cc.f DIS/nc_dis.f DIS/xcoeffunc_cc.f DIS/xpdfs.f DIS/dis_xsec.f DIS/FFCoefficientFunctions.f DIS/integrands_nc.f DIS/PDFparam.f DIS/xcoeffunc_nc.f DIS/ZMCoefficientFunctions.f

EXTRA_DIST = commons

HELL = HELL/src/hell.cc HELL/src/DGLAPevol.cc HELL/src/hell_wrapper.cc

libAPFELCore_la_SOURCES = $(INIT) $(SETTINGS) $(INTER) $(SPLITTING_FUNCTIONS) $(COUPLINGS) $(EVOLUTIONQCD) $(EVOLUTIONQED) $(PDF) $(UTILITIES) $(DIS) $(HELL)

AM_CXXFLAGS += -I$(shell lhapdf-config --incdir)
