.PHONY: skala-reconstruction-test skala-primitives-test skala-partition-probe skala-periodic-test skala-lebedev-test
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

skala-periodic-test: libpaw.a
	mkdir -p unit-tests
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/periodic_partition.x $(BASEDIR)/tests/unittests/skala_reconstruction/periodic_partition.f90 libpaw.a $(LIBS)
	./unit-tests/periodic_partition.x

skala-lebedev-test: libpaw.a
	mkdir -p unit-tests
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/lebedev_exactness.x $(BASEDIR)/tests/unittests/skala_reconstruction/lebedev_exactness.f90 libpaw.a $(LIBS)
	./unit-tests/lebedev_exactness.x
	@if ./unit-tests/lebedev_exactness.x 66 > unit-tests/lebedev_invalid.log 2>&1; then exit 1; fi
	grep -q 'Lebedev minimum exactness must lie in 1..65' unit-tests/lebedev_invalid.log
