# CC = gcc
CC = icx
CFLAGS = -O3 -I ../
# FORT = gfortran
FORT = mpiifx
FFLAGS = -O3 -I ../../fortran/
LFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lifcore

all: test

test: test.out
	./test.out

test.out: check.o ../libcpfapack.a ../../fortran/libpfapack.a
	$(CC) $(CFLAGS) test.c -o test.out check.o ../libcpfapack.a \
	    ../../fortran/libpfapack.a $(LFLAGS)

clean:
	rm -f *.o
	rm -f *.out

.SUFFIXES: .o .f90

.f90.o:
	$(FORT) $(FFLAGS) -c $< -o $@
