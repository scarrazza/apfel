AM_FCFLAGS  = $(FCLIBS)
AM_LDFLAGS = $(top_builddir)/lib/libAPFEL.la -L$(shell lhapdf-config --libdir) -lLHAPDF

noinst_PROGRAMS = Tabulation SumRules LHgridProduction LHgridDerivativeProduction Luminosity DISObservables TabulationStep TabulationFrag TabulationExternal ListFunctions
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

noinst_PROGRAMS += TabulationCxx SumRulesCxx LHgridProductionCxx LHgridDerivativeProductionCxx LuminosityCxx Timing DISObservablesCxx TabulationStepCxx TabulationFragCxx ListFunctionsCxx
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

EXTRA_DIST = TabulationPython SumRulesPython LHgridProductionPython LuminosityPython
