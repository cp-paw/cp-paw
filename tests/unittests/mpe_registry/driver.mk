.PHONY: mpe-registry-test
mpe-registry-test: paw.x
	mkdir -p unit-tests
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/mpe_registry.x $(BASEDIR)/tests/unittests/mpe_registry/registry.f90 $(ARGLIST) $(LIBS)
	$(MPIEXEC) -np 3 ./unit-tests/mpe_registry.x
