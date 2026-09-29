.PHONY: skala-reconstruction-test skala-primitives-test skala-partition-probe
skala-reconstruction-test: libpaw.a
	mkdir -p unit-tests
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/skala_reconstruction.x $(BASEDIR)/tests/unittests/skala_reconstruction/skala_reconstruction.f90 libpaw.a $(LIBS)
	./unit-tests/skala_reconstruction.x

skala-primitives-test: libpaw.a
	mkdir -p unit-tests
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/skala_primitives.x $(BASEDIR)/tests/unittests/skala_primitives/skala_primitives.f90 libpaw.a $(LIBS)
	./unit-tests/skala_primitives.x

skala-partition-probe: libpaw.a
	mkdir -p unit-tests
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/partition_measure.x $(BASEDIR)/tests/unittests/skala_reconstruction/partition_measure.f90 libpaw.a $(LIBS)
