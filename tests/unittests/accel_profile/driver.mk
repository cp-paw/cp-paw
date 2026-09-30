.PHONY: accel-profile-report-test
accel-profile-report-test: libpaw.a
	mkdir -p unit-tests
	$(F90PP) $(CPPFLAGS) < $(BASEDIR)/tests/unittests/accel_profile/report.f90 > unit-tests/profile_report.f90
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/profile_report.x unit-tests/profile_report.f90 libpaw.a $(LIBS)
	python3 $(BASEDIR)/tests/unittests/accel_profile/check_report.py unit-tests/profile_report.x
