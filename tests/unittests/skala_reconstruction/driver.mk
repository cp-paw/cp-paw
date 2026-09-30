.PHONY: skala-reconstruction-test skala-primitives-test skala-partition-probe skala-periodic-test skala-lebedev-test skala-partition-cache-test skala-interpolation-test

skala-interpolation-test: libpaw.a
	mkdir -p unit-tests
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/interpolation.x $(BASEDIR)/tests/unittests/skala_reconstruction/interpolation.f90 libpaw.a $(LIBS)
	./unit-tests/interpolation.x

skala-reconstruction-test: libpaw.a
	mkdir -p unit-tests
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/skala_reconstruction.x $(BASEDIR)/tests/unittests/skala_reconstruction/skala_reconstruction.f90 libpaw.a $(LIBS)
	./unit-tests/skala_reconstruction.x
	CPPAW_SKALA_SOURCE_CACHE_MB=0 ./unit-tests/skala_reconstruction.x limit > unit-tests/source_cache_off.log
	grep -qx '0' unit-tests/source_cache_off.log
	CPPAW_SKALA_SOURCE_CACHE_MB=256 ./unit-tests/skala_reconstruction.x limit > unit-tests/source_cache_limit.log
	grep -qx '268435456' unit-tests/source_cache_limit.log
	@if GFORTRAN_UNBUFFERED_ALL=1 CPPAW_SKALA_SOURCE_CACHE_MB=invalid ./unit-tests/skala_reconstruction.x limit > unit-tests/source_cache_invalid.log 2>&1; then exit 1; fi
	grep -q 'CPPAW_SKALA_SOURCE_CACHE_MB MUST BE A NONNEGATIVE INTEGER' unit-tests/source_cache_invalid.log

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

skala-partition-cache-test: libpaw.a
	mkdir -p unit-tests
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/partition_cache.x $(BASEDIR)/tests/unittests/skala_reconstruction/partition_cache.f90 libpaw.a $(LIBS)
	./unit-tests/partition_cache.x
	CPPAW_SKALA_PARTITION_CACHE_MB=0 ./unit-tests/partition_cache.x limit > unit-tests/partition_cache_off.log
	grep -qx '0' unit-tests/partition_cache_off.log
	CPPAW_SKALA_PARTITION_CACHE_MB=256 ./unit-tests/partition_cache.x limit > unit-tests/partition_cache_limit.log
	grep -qx '268435456' unit-tests/partition_cache_limit.log
	@if GFORTRAN_UNBUFFERED_ALL=1 CPPAW_SKALA_PARTITION_CACHE_MB=invalid ./unit-tests/partition_cache.x limit > unit-tests/partition_cache_invalid.log 2>&1; then exit 1; fi
	grep -q 'CPPAW_SKALA_PARTITION_CACHE_MB MUST BE A NONNEGATIVE INTEGER' unit-tests/partition_cache_invalid.log

skala-lebedev-test: libpaw.a
	mkdir -p unit-tests
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/lebedev_exactness.x $(BASEDIR)/tests/unittests/skala_reconstruction/lebedev_exactness.f90 libpaw.a $(LIBS)
	./unit-tests/lebedev_exactness.x
	@if GFORTRAN_UNBUFFERED_ALL=1 ./unit-tests/lebedev_exactness.x 66 > unit-tests/lebedev_invalid.log 2>&1; then exit 1; fi
	grep -q 'Lebedev minimum exactness must lie in 1..65' unit-tests/lebedev_invalid.log
