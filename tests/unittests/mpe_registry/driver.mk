MPIEXEC ?= mpirun
.PHONY: mpe-registry-test
mpe-registry-test: paw.x
	mkdir -p unit-tests
	$(F90PP) $(CPPFLAGS) < $(BASEDIR)/tests/unittests/mpe_registry/registry.f90 > unit-tests/mpe_registry.f90
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/mpe_registry.x unit-tests/mpe_registry.f90 $(ARGLIST) $(LIBS)
	$(MPIEXEC) -np 3 ./unit-tests/mpe_registry.x

.PHONY: skala-point-exchange-test
skala-point-exchange-test: paw.x
	mkdir -p unit-tests
	$(F90PP) $(CPPFLAGS) < $(BASEDIR)/tests/unittests/mpe_registry/points.f90 > unit-tests/skala_points.f90
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/skala_points.x unit-tests/skala_points.f90 $(ARGLIST) $(LIBS)
	$(MPIEXEC) -np 3 ./unit-tests/skala_points.x
