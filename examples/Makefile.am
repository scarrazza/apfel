AM_FCFLAGS  = $(FCLIBS)
AM_LDFLAGS = $(top_builddir)/lib/libAPFEL.la $(LHAPDF_LDFLAGS)

noinst_PROGRAMS = Tabulation SumRules LHgridProduction LHgridDerivativeProduction Luminosity DISObservables TabulationStep TabulationFrag TabulationExternal ListFunctions CheckAPFEL
Tabulation_SOURCES = Tabulation.f
SumRules_SOURCES = SumRules.f
LHgridProduction_SOURCES = LHgridProduction.f
LHgridDerivativeProduction_SOURCES = LHgridDerivativeProduction.f
Luminosity_SOURCES = Luminosity.f
DISObservables_SOURCES = DISObservables.f
TabulationStep_SOURCES = TabulationStep.f
TabulationFrag_SOURCES = TabulationFrag.f
TabulationExternal_SOURCES = TabulationExternal.f
ListFunctions_SOURCES = ListFunctions.f
CheckAPFEL_SOURCES = CheckAPFEL.f

noinst_PROGRAMS += TabulationCxx SumRulesCxx LHgridProductionCxx LHgridDerivativeProductionCxx LuminosityCxx Timing DISObservablesCxx TabulationStepCxx TabulationFragCxx ListFunctionsCxx CheckAPFELCxx
TabulationCxx_SOURCES = TabulationCxx.cc
SumRulesCxx_SOURCES = SumRulesCxx.cc
LHgridProductionCxx_SOURCES = LHgridProductionCxx.cc
LHgridDerivativeProductionCxx_SOURCES = LHgridDerivativeProductionCxx.cc
Timing_SOURCES = Timing.cc
LuminosityCxx_SOURCES = LuminosityCxx.cc
DISObservablesCxx_SOURCES = DISObservablesCxx.cc
TabulationStepCxx_SOURCES = TabulationStepCxx.cc
TabulationFragCxx_SOURCES = TabulationFragCxx.cc
ListFunctionsCxx_SOURCES = ListFunctionsCxx.cc
CheckAPFELCxx_SOURCES = CheckAPFELCxx.cc

EXTRA_DIST = TabulationPython SumRulesPython LHgridProductionPython LuminosityPython
