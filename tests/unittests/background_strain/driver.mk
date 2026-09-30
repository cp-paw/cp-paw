.PHONY: background-strain-test
background-strain-test: paw.x
	mkdir -p unit-tests
	$(F90PP) $(CPPFLAGS) < $(BASEDIR)/tests/unittests/background_strain/background_strain.f90 > unit-tests/background_strain.f90
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/background_strain.x unit-tests/background_strain.f90 $(ARGLIST) $(LIBS)
	./unit-tests/background_strain.x
