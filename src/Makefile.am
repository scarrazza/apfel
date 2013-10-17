noinst_LTLIBRARIES = libAPFELCore.la

INIT = 	apfelversion.cc WelcomeMessage.f initParameters.f initGrid.f initIntegralsQCD.f \
	initIntegralsQED.f initIntegralsMatching.f initPDFs.f \
	InitializeAPFEL.f EvolveAPFEL.f

SETTINGS = SetAlphaQCDRef.f SetAlphaQEDRef.f SetQLimits.f \
	SetPerturbativeOrder.f SetFFNS.f SetGridParameters.f \
	SetMSbarMasses.f SetNumberOfGrids.f SetPDFSet.f SetPoleMasses.f \
	SetRenFacRatio.f SetTheory.f SetVFNS.f SetMaxFlavourPDFs.f \
	SetMaxFlavourAlpha.f SetReplica.f AlphaQCD.f AlphaQED.f LHAPDFgrid.f

INTER = interpolants.f

SPLITTING_FUNCTIONS = SplittingFunctions.f

COUPLINGS = a_QCD.f a_QED.f lambda.f

EVOLUTIONQCD = EvolutionOperatorsQCD.f IdentityOperatorsQCD.f \
	EqualOperatorsQCDnf.f integrandsQCD.f integralsQCD.f RSLintegralsQCD.f \
	odeintnsQCD.f odeintsgQCD.f MatchingConditions.f integrandsMatching.f \
	integralsMatching.f RSLintegralsMatching.f MatchPDFs.f

EVOLUTIONQED = EvolutionOperatorsQED.f IdentityOperatorsQED.f \
	EqualOperatorsQEDnf.f integrandsQED.f integralsQED.f RSLintegralsQED.f \
	odeintnsQED.f odeintsgQED.f

PDF = PDFphys2evQCD.f PDFevQCD2phys.f PDFphys2evQED.f PDFevQED2phys.f \
	EvolvePDFs.f switchGluonPhoton.f EvolveQCD.f EvolveQED.f toyLHPDFs.f \
	xPDF.f xgamma.f NPDF.f Ngamma.f LUMI.f

UTILITIES = dgauss.f wgplg.f ddilog.f

libAPFELCore_la_SOURCES = $(INIT) $(SETTINGS) $(INTER) $(SPLITTING_FUNCTIONS) \
	$(COUPLINGS) $(EVOLUTIONQCD) $(EVOLUTIONQED) $(PDF) $(UTILITIES)
