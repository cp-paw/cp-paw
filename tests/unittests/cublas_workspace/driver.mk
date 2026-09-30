.PHONY: cublas-workspace-build cublas-workspace-test

cublas-workspace-build: libpaw.a
	mkdir -p unit-tests
	$(F90PP) $(CPPFLAGS) < $(BASEDIR)/tests/unittests/cublas_workspace/cublas_workspace.f90 > unit-tests/cublas_workspace.f90
	$(FC) $(FCFLAGS) -I. -c unit-tests/cublas_workspace.f90 -o unit-tests/cublas_workspace.o
	$(LD) $(LDFLAGS) -o unit-tests/cublas_workspace.x unit-tests/cublas_workspace.o libpaw.a $(LIBS)

cublas-workspace-test: cublas-workspace-build
	CPPAW_CUBLAS_ACC=1 CPPAW_CUBLAS_ACC_MINFLOP=1 CPPAW_CUBLAS_ACC_MATMUL_MINFLOP=1 CPPAW_CUBLAS_ACC_ADDPRODUCT_MINFLOP=1 ./unit-tests/cublas_workspace.x
