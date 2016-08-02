EXTRA_DIST = apfel.i

## Bundling SWIG-built sources with the tarball for SL4 compat
EXTRA_DIST += apfel_wrap.cc apfel.py
apfel_wrap.cc apfel.py: apfel.i
	swig -c++ -python -I$(HEPMCINCPATH) -I$(top_srcdir)/include -o apfel_wrap.cc $<

AM_CXXFLAGS += -Wno-long-long

all-local: apfel_wrap.cc
	$(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix)

uninstall-local:
	rm -rf $(DESTDIR)$(libdir)/python*/site-packages/apfel*

clean-local:
	$(PYTHON) setup.py clean --all
	@rm -f $(top_builddir)/*.pyc

distclean-local:
	@rm -rf apfel.egg-info
	@rm -rf dist
	@rm -f apfel.py
	@rm -f apfel_wrap.cc

