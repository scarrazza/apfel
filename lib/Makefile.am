lib_LTLIBRARIES = libAPFEL.la

libAPFEL_la_SOURCES = 

libAPFEL_la_LDFLAGS = $(VERSIONFLAGS)

libAPFEL_la_LIBADD = $(FCLIBS) $(top_builddir)/src/libAPFELCore.la \
	$(top_builddir)/ccwrap/libAPFELCCWrapEvol.la \
	$(top_builddir)/ccwrap/libAPFELCCWrapObs.la 

## Hack to get linking done by CXX
libAPFEL_la_LINK = $(LIBTOOL) --tag=CXX --mode=link $(CXX) $(AM_CXXFLAGS) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@
