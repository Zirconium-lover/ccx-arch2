all: ccx_2.22

CFLAGS = -Wall -O3 -march=native -mtune=native -DARCH="Linux" -DPARDISO -DARPACK -DMATRIXSTORAGE -DNETWORKOUT -fopenmp -I/opt/intel/oneapi/mkl/latest/include
FFLAGS = -Wall -O3 -march=native -mtune=native -fopenmp

CC = gcc
FC = gfortran

.c.o:
	$(CC) $(CFLAGS) -c $<
.f.o:
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.22.c
OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)

ccx_2.22.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

ccx_2.22: ccx_2.22.a ccx_2.22.c
	$(CC) $(CFLAGS) -c ccx_2.22.c
	$(FC) $(FFLAGS) -o $@ ccx_2.22.o ccx_2.22.a \
	/mnt/pool/6/railyasov/modelingTensor/tools/arpack-ng/build/libarpack.a \
	-L/opt/intel/oneapi/mkl/latest/lib/intel64 \
	-lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -fopenmp -lpthread -lm -ldl

clean:
	rm -f *.o *.mod *.a ccx_2.22 ccx_2.22.o
