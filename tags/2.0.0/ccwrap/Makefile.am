EXTRA_DIST = Doxyfile

## C++ wrapper library to be combined into libAPFEL
noinst_LTLIBRARIES = libAPFELCCWrap.la
libAPFELCCWrap_la_SOURCES = APFELfw.f APFEL.cc

## Deal with the Doxygen stuff
if WITH_DOXYGEN
.PHONY = dox
dox:
	$(DOXYGEN) $(top_srcdir)/ccwrap/Doxyfile || true
all-local: dox
	@true
install-data-local:
	$(install_sh_DATA) doxy/html/* $(DESTDIR)$(pkgdatadir)/doc/html
mostlyclean-local:
	rm -rf $(top_builddir)/ccwrap/doxy
uninstall-local:
	rm -rf $(DESTDIR)$(pkgdatadir)/doc; \
	test ! -d $(DESTDIR)$(pkgdatadir) || rmdir --ignore-fail-on-non-empty $(DESTDIR)$(pkgdatadir)
endif
