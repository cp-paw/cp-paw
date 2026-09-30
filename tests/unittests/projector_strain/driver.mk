.PHONY: projector-strain-test
projector-strain-test: paw.x
	mkdir -p unit-tests
	$(F90PP) $(CPPFLAGS) < $(BASEDIR)/tests/unittests/projector_strain/projector_strain.f90 > unit-tests/projector_strain.f90
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/projector_strain.x unit-tests/projector_strain.f90 $(ARGLIST) $(LIBS)
	./unit-tests/projector_strain.x
