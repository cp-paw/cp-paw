.PHONY: setup-reciprocal-test
setup-reciprocal-test: paw.x
	mkdir -p unit-tests
	$(F90PP) $(CPPFLAGS) < $(BASEDIR)/tests/unittests/setup_reciprocal/setup_reciprocal.f90 > unit-tests/setup_reciprocal.f90
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/setup_reciprocal.x unit-tests/setup_reciprocal.f90 $(ARGLIST) $(LIBS)
	./unit-tests/setup_reciprocal.x
