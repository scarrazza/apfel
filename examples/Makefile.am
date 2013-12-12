AM_FCFLAGS  = $(FCLIBS)
AM_LDFLAGS = $(top_builddir)/lib/libAPFEL.la -L$(shell lhapdf-config --libdir) -lLHAPDF

noinst_PROGRAMS = Tabulation SumRules LHgridProduction Luminosity DISObservables TabulationStep
Tabulation_SOURCES = Tabulation.f
SumRules_SOURCES = SumRules.f
LHgridProduction_SOURCES = LHgridProduction.f
Luminosity_SOURCES = Luminosity.f
DISObservables_SOURCES = DISObservables.f
TabulationStep_SOURCES = TabulationStep.f

noinst_PROGRAMS += TabulationCxx SumRulesCxx LHgridProductionCxx LuminosityCxx Timing DISObservablesCxx
TabulationCxx_SOURCES = TabulationCxx.cc
SumRulesCxx_SOURCES = SumRulesCxx.cc
LHgridProductionCxx_SOURCES = LHgridProductionCxx.cc
Timing_SOURCES = Timing.cc
LuminosityCxx_SOURCES = LuminosityCxx.cc
DISObservablesCxx_SOURCES = DISObservablesCxx.cc

EXTRA_DIST = TabulationPython SumRulesPython LHgridProductionPython LuminosityPython
