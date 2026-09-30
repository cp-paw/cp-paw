.PHONY: skala-reconstruction-test skala-source-tiles-test skala-primitives-test skala-partition-probe skala-periodic-test skala-lebedev-test skala-partition-cache-test skala-interpolation-test

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

# Requires an OpenACC build and a visible NVIDIA GPU. The row oracle covers
# spin, existing matrix contents, partial caches and non-full final tiles.
skala-source-tiles-test: export OMP_NUM_THREADS := 1
skala-source-tiles-test: export OPENBLAS_NUM_THREADS := 1
skala-source-tiles-test: skala-reconstruction-test
	@for rows in 1 7 128 512 2048; do \
	  CPPAW_SKALA_SOURCE_BACK_ACC=1 CPPAW_SKALA_SOURCE_BACK_ACC_MIN_ROWS=1 \
	  CPPAW_SKALA_SOURCE_BACK_ACC_MB=256 CPPAW_SKALA_TEST_REQUIRE_SOURCE_ACC=1 \
	  CPPAW_SKALA_SOURCE_BACK_ACC_TILE_ROWS=$$rows ./unit-tests/skala_reconstruction.x \
	    > unit-tests/source_tile_$$rows.log 2>&1 || exit 1; \
	done
	CPPAW_SKALA_SOURCE_BACK_ACC=1 CPPAW_SKALA_SOURCE_BACK_ACC_MB=0 CPPAW_SKALA_TEST_REQUIRE_SOURCE_ACC=0 CPPAW_SKALA_SOURCE_BACK_ACC_TILE_ROWS=2048 ./unit-tests/skala_reconstruction.x > unit-tests/source_tile_zero_budget.log 2>&1
	grep -qx 'SOURCE GPU BATCH ROWS 0' unit-tests/source_tile_zero_budget.log
	@for rows in 0 -1 invalid 999999999999 ''; do \
	  if CPPAW_SKALA_SOURCE_BACK_ACC=1 CPPAW_SKALA_SOURCE_BACK_ACC_TILE_ROWS="$$rows" \
	    ./unit-tests/skala_reconstruction.x > unit-tests/source_tile_invalid.log 2>&1; then exit 1; fi; \
	  grep -q 'CPPAW_SKALA_SOURCE_BACK_ACC_TILE_ROWS MUST BE A POSITIVE INTEGER' \
	    unit-tests/source_tile_invalid.log || exit 1; \
	done

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
