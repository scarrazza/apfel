CC = g++
#FF = gfortran

PKG = hell

CFLAGS = -O2 -Wall -I.

obj = objects

OBJECTS = $(obj)/hell.o $(obj)/DGLAPevol.o $(obj)/gammaNLO.o



all: test_hell test_evolution
lib: lib$(PKG).a


test%: test%.cc lib$(PKG).a
	$(CC) $(CFLAGS) -o $@ $< lib$(PKG).a -lgsl -lgslcblas

lib$(PKG).a: $(OBJECTS)
	ar -r lib$(PKG).a $^

$(obj)/%.o: src/%.cc include/%.hh include/math/matrix.hh
	$(CC) $(CFLAGS) -o $@ -c $<

#$(obj)/%.o: math/%.cc math/%.hh
#	$(CC) $(CFLAGS) -o $@ -c $<




clean: clean-tilde
	rm -f test_hell
	rm -f test_evolution
	rm -f lib$(PKG).a
	rm -f $(obj)/*.o

clean-tilde:
	rm -f *~
	rm -f src/*~
	rm -f src/math/*~
	rm -f include/*~
	rm -f include/math/*~


.PHONY: clean-tilde clean lib