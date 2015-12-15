lib_LTLIBRARIES = libAPFELevol.la
noinst_LTLIBRARIES = libAPFELCore.la

INIT = 	Evolution/apfelversion.cc Evolution/lhapdfwrapper.cc Evolution/WelcomeMessage.f Evolution/initParameters.f \
	Evolution/initGrid.f Evolution/initGridAlpha.f Evolution/initIntegralsQCD.f \
	Evolution/initIntegralsQCDRes.f Evolution/initIntegralsQED.f Evolution/initIntegralsMatching.f \
	Evolution/initPDFs.f Evolution/InitializeAPFEL.f Evolution/EvolveAPFEL.f Evolution/DeriveAPFEL.f \
	Evolution/initHELL.f Evolution/TruncatedEvolveAPFEL.f Evolution/ExponentiatedEvolveAPFEL.f \
	Evolution/ReportParameters.f Evolution/CheckAPFEL.f Evolution/ComputeHeavyQuarkThresholds.f \
	Evolution/CachePDFsAPFEL.f

SETTINGS = Evolution/SetAlphaQCDRef.f Evolution/SetAlphaQEDRef.f Evolution/SetQLimits.f \
	Evolution/SetPerturbativeOrder.f Evolution/SetFFNS.f Evolution/SetGridParameters.f Evolution/SetExternalGrid.f \
	Evolution/SetMSbarMasses.f Evolution/SetNumberOfGrids.f Evolution/SetPDFSet.f Evolution/SetPoleMasses.f \
	Evolution/SetRenFacRatio.f Evolution/SetTheory.f Evolution/SetVFNS.f Evolution/SetMaxFlavourPDFs.f \
	Evolution/SetMaxFlavourAlpha.f Evolution/SetReplica.f Evolution/AlphaQCD.f Evolution/AlphaQED.f Evolution/HeavyQuarkMass.f \
	Evolution/LHAPDFgrid.f Evolution/CleanUp.f Evolution/EnableWelcomeMessage.f Evolution/EnableEvolutionOperator.f \
	Evolution/ExternalEvolutionOperator.f Evolution/LockGrids.f Evolution/SetTimeLikeEvolution.f Evolution/SetAlphaEvolution.f \
	Evolution/SetLambdaQCDRef.f Evolution/SetPDFEvolution.f Evolution/LHAPDFgridDerivative.f \
	Evolution/SetSmallxResummation.f Evolution/SetFastEvolution.f Evolution/EnableMassRunning.f Evolution/ListFunctions.f \
	Evolution/GetPerturbativeOrder.f Evolution/SetTauMass.f Evolution/EnableLeptonEvolution.f Evolution/SetLHgridParameters.f \
	Evolution/SetEpsilonTruncation.f Evolution/GetMuF.f Evolution/SetPolarizedEvolution.f Evolution/SetMassMatchingScales.f \
	Evolution/SetMassScaleReference.f

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
	Evolution/xGrid.f Evolution/nIntervals.f Evolution/xLepton.f Evolution/xPDFxQ.f

UTILITIES = core/dgauss.f core/wgplg.f core/ddilog.f core/gamma.f core/beta.f core/hplog.f

DIS =   DIS/InitializeAPFEL_DIS.f DIS/initParametersDIS.f DIS/SetMassScheme.f DIS/initIntegralsDIS.f \
	DIS/ZeroMassCoefficientFunctions.f DIS/RSLintegralsDIS.f DIS/ComputeStructureFunctionsAPFEL.f DIS/F2total.f \
	DIS/F2light.f DIS/F2charm.f DIS/F2bottom.f DIS/F2top.f DIS/FLtotal.f DIS/FLlight.f DIS/FLcharm.f \
	DIS/FLbottom.f DIS/FLtop.f DIS/F3total.f DIS/F3light.f DIS/F3charm.f DIS/F3bottom.f DIS/F3top.f \
	DIS/integrandsDIS.f DIS/SetProcessDIS.f DIS/ComputeChargesDIS.f DIS/SetPolarizationDIS.f DIS/SetProjectileDIS.f \
	DIS/SetTargetDIS.f DIS/MassiveCoefficientFunctions.f DIS/MassiveZeroCoefficientFunctions.f DIS/ComputeDISOperators.f \
	DIS/ConvolutePDFsWithDISOperators.f DIS/JoinDISOperators.f DIS/ConvoluteEvolutionWithDISOperators.f \
	DIS/ExternalDISOperator.f DIS/SetCKM.f DIS/SetGFermi.f DIS/SetProtonMass.f DIS/SetSin2ThetaW.f \
	DIS/SetWMass.f DIS/SetZMass.f DIS/FKSimulator.f DIS/SetFKObservable.f DIS/GetFKObservable.f DIS/FKObservables.f \
	DIS/EnableTargetMassCorrections.f DIS/GetCKM.f DIS/GetGFermi.f DIS/GetProtonMass.f DIS/GetSin2ThetaW.f \
	DIS/GetWMass.f DIS/GetZMass.f DIS/RSLintegralsSIA.f DIS/integrandsSIA.f DIS/initIntegralsSIA.f DIS/SelectCharge.f \
	DIS/SetRenQRatio.f DIS/SetFacQRatio.f DIS/GetSIATotalCrossSection.f DIS/IncludeScaleVariation.f \
	DIS/ExcludeScaleVariation.f DIS/ReportParametersDIS.f DIS/EnableDampingFONLL.f DIS/SetPropagatorCorrection.f \
	DIS/SetEWCouplings.f DIS/EnableDynamicalScaleVariations.f DIS/EnableIntrinsicCharm.f DIS/IncludeIntrinsicCharm.f

FTDY = FTDY/src/ComputeFKTables.f FTDY/src/ComputeHardCrossSectionsDY.f FTDY/src/gauleg.f FTDY/src/initDYcouplings.f FTDY/src/initxGridDY.f \
	FTDY/src/readcDY.f FTDY/src/ReadDataFile.f FTDY/src/sigmafk_dy.f FTDY/src/wcDYlo.f FTDY/src/wcDYqg.f FTDY/src/wcDYqq.f FTDY/src/wcoeffDY.f \
	FTDY/src/writecDY.f FTDY/src/writeFK.f FTDY/src/xlineigen.f FTDY/src/PDFphys2ev.f

EXTRA_DIST = commons FTDY/commons

HELL = HELL/src/hell.cc HELL/src/hell_wrapper.cc

libAPFELCore_la_SOURCES = $(INIT) $(SETTINGS) $(INTER) $(SPLITTING_FUNCTIONS) $(COUPLINGS) $(EVOLUTIONQCD) $(EVOLUTIONQED) $(PDF) $(UTILITIES) $(DIS) $(HELL) $(FTDY)

libAPFELevol_la_LIBADD = $(FCLIBS) $(top_builddir)/ccwrap/libAPFELCCWrapEvol.la
libAPFELevol_la_SOURCES = $(INIT) $(SETTINGS) $(INTER) $(SPLITTING_FUNCTIONS) $(COUPLINGS) $(EVOLUTIONQCD) $(EVOLUTIONQED) $(PDF) $(UTILITIES) $(HELL)

AM_CXXFLAGS += $(LHAPDF_CPPFLAGS) $(LHAPDF_CXXFLAGS)
