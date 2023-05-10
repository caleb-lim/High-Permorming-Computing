# here the make accepts the compiler, gpu, calctype
# will add a make print to indicate what types of compilations
# are enabled.
COMPILERTYPE ?= GCC
GPUTYPE ?= NVIDIA
OPTLEVEL ?= 2
VISUALTYPE ?= STANDARD
PROFILER ?= OFF

# optmisation flags
OPTFLAGS = -O$(OPTLEVEL)
# profiling flags if desired 
ifeq ($(PROFILER), ON)
	OPTFLAGS += -pg -g
endif

# formatting characters for info output
NULL := 
TAB := $(NULL)  $(NULL)

# lets define a bunch of compilers
GCC = gcc
GCCCXX = g++
GCCFORT = gfortran
CLANG = clang
CLANGCXX = clang++
CLANGFORT = flang
AOMP = aompcc
AOMPCXX = aompcc
CRAYCC = cc 
CRAYCXX = CC
CRAYFORT = ftn 
INTELCC = icc 
INTELCXX = icpc
INTELFORT = ifort

# C flags. Currently only really need to enforce c11 standard
GCCCFLAGS = -std=c11

# fortran flags 
GCCFFLAGS = -cpp -ffixed-line-length-none -dM
INTELFFLAGS = -cpp -extend-source -D_INTELFTN

# cuda
NCC = nvcc
NCXX = nvcc
CUDA_FLAGS = -DUSECUDA

# hip
HCC = hipcc
HCXX = hipcc
HIP_FLAGS = -DUSEHIP
OCL_FLAGS = -lOpenCL -DUSEOPENCL

# mpi
MPICC = mpicc
MPICXX = mpicxx
MPIFORT = mpif90

# openmp 
GCCOMP_FLAGS = -fopenmp -DUSEOPENMP
GCCOMPTARGET_FLAGS = -fopenmp -DUSEOPENMPTARGET
INTELOMP_FLAGS = -qopenmp -DUSEOPENMP
INTELMPTARGET_FLAGS = -qopenmp -DUSEOPENMPTARGET
GCCOACC_FLAGS = -fopenacc -fopt-info-optimized-omp -DUSEOPENACC
INTELOACC_FLAGS = -qopenacc -fopt-info-optimized-omp -DUSEOPENACC

#now set the default compilers
CC = $(GCC)
CXX = $(GCCCXX)
FORT = $(GCCFORT)
GPUCC = $(NCC)
GPUCXX = $(NCXX)
FFLAGS = $(GCCFFLAGS)
OMP_FLAGS = $(GCCOMP_FLAGS)
OMPTARGET_FLAGS = $(GCCOMPTARGET_FLAGS)
OACC_FLAGS = $(GCCOACC_FLAGS)
CFLAGS = $(GCCCFLAGS)

# change if required
ifeq ($(COMPILERTYPE), CRAY-GNU)
	CC = $(CRAYCC)
	CXX = $(CRAYCXX)
	FORT = $(CRAYFORT)
endif
ifeq ($(COMPILERTYPE), CLANG)
	CC = $(CLANG)
	CXX = $(CLANGCXX)
	FORT = $(CLANGFORT)        
endif
ifeq ($(COMPILERTYPE), AOMP)
	CC = $(AOMP)
	CXX = $(AOMPCXX)
endif
ifeq ($(COMPILERTYPE), CRAY)
	CC = $(CRAYCC)
	CXX = $(CRAYCXX)
	CFLAGS =
	FFLAGS = -eZ -ffree
endif
ifeq ($(COMPILERTYPE), INTEL)
	CC = $(INTELCC)
	CXX = $(INTELCXX)
	FORT = $(INTELFORT)
	FFLAGS = $(INTELFFLAGS)
	OMP_FLAGS = $(INTELOMP_FLAGS)
	OMPTARGET_FLAGS = $(INTELOMPTARGET_FLAGS)
endif
ifeq ($(COMPILERTYPE), CRAY-INTEL)
	CC = $(CRAYCC)
	CXX = $(CRAYCXX)
	FORT = $(CRAYFORT)
	FFLAGS = $(INTELFFLAGS)
	OMP_FLAGS = $(INTELOMP_FLAGS)
	OMPTARGET_FLAGS = $(INTELOMPTARGET_FLAGS)
endif
ifeq ($(GPUTYPE), AMD)
	GPUCC = $(HCC)
	GPUCCXX = $(HCXX)
endif

# compiler used for openmp
OMPCC = $(CC)
OMPCXX = $(CXX)
OMPFORT = $(FORT)

# compilers used for openacc
OACCCC = $(CC)
OACCCXX = $(CXX)

# compilers used for opencl
OCLC = $(CCGPU)
OCLCXX = $(CXXGPU)

VISUALFLAGS =
# here this flag may need to be updated with PNG location
PNGDIR ?= /opt/homebrew/
PNGLIB ?= -lpng
ifeq ($(VISUALTYPE), PNG)
	VISUALFLAGS += $(PNGLIB) -I$(PNGDIR)/include/ -L$(PNGDIR)/lib/ -DUSEPNG
endif

COMMONFLAGS = $(OPTFLAGS) $(VISUALFLAGS)

.PHONY : dirs allinfo configinfo buildinfo makecommands clean
.PHONY : cpu_serial cpu_serial_cc cpu_serial_fort 
.PHONY : cpu_serial_expanded cpu_serial_expanded_cc cpu_serial_expanded_fort
.PHONY : cpu_openmp cpu_openmp_loop cpu_openmp_task 
.PHONY : cpu_openmp_cc cpu_openmp_loop_cc cpu_openmp_task_cc 
.PHONY : cpu_openmp_fort cpu_openmp_loop_fort cpu_openmp_task_fort
.PHONY : cpu_mpi cpu_mpi_cc cpu_mpi_fort
.PHONY : gpu_openmp gpu_openacc hip_monogpu

all : dirs buildinfo cpu_serial cpu_serial_expanded 

allinfo : configinfo makecommands buildinfo 

# have it spit out information about current configuration
configinfo : 
	$(info )
	$(info ========================================)
	$(info Compiler options:)
	$(info ========================================)
	$(info Allowed COMPILERTYPE: GCC, CLANG, INTEL, AOMP, CRAY, CRAY-GNU, CRAY-INTEL)
	$(info Allowed GPUTYPE: NVIDIA, AMD)
	$(info )
	$(info Other options:)
	$(info ----------------------------------------)
	$(info One can also specify optmisation level with OPTLEVEL=)
	$(info one can use 0, 1, 2, 3, or fast)
	$(info )
	$(info One can also turn on profiling flags using PROFILER=ON)
	$(info )
	$(info One can also specify using PNG visualisation of GOL with VISUALTYPE=PNG)
	$(info but does require setting the PNGDIR to point to the base directory)
	$(info where include/png.h and lib/libpng are located.)
	$(info if the png library has a non-standard name like /usr/lib/libpng16.so.16)
	$(info set PNGLIB= to the non-standard name)
	$(info ========================================)
	$(info )


# list current build info given the commands make was passed
buildinfo :
	$(info )
	$(info ========================================)
	$(info Current compilers to be used:)
	$(info ========================================)
	$(info Compiling with $(CC) $(CXX) $(FORT) for CPU focused codes)
	$(info Compiling with $(MPICC) $(MPICXX) $(MPIFORT) for MPI-CPU focused codes)
	$(info Compiling with $(GPUCC) $(GPUCXX) for GPU focused codes)
	$(info Compiling with $(OMPCC) $(OMPCXX) $(OMPFORT) for OpenMP directive GPU focused codes)
	$(info Compiling with $(OACCCC) $(OACCCXX) $(OACCFORT) for OpenACC directive GPU focused codes)
	$(info ========================================)
	$(info )

# list current build info given the commands make was passed
makecommands :
	$(info )
	$(info ========================================)
	$(info Make commands:)
	$(info ========================================)
	$(info Make is configured so that the following can be compiled if provided this argument)
	$(info )
	$(info cpu_serial : compiles both C and Fortran verisons of the standard serial code. )
	$(info $(TAB)Sources : 01_gol_cpu_serial.c or 01_gol_cpu_serial_fort.f90)
	$(info cpu_serial_expanded : compiles C verison of serial code setup to allow for easy expansion of rules. )
	$(info $(TAB)Sources : 01_gol_cpu_serial_expanded.c)
	$(info )
	$(info cpu_openmp : compiles both C and Fortran verisons cpu parallelised openmp codes.)
	$(info Will try compiling both loop and task based parallelism )
	$(info $(TAB)Sources : 02_gol_cpu_openmp_loop.c or 02_gol_cpu_openmp_loop_fort.f90)
	$(info $(TAB)          02_gol_cpu_openmp_task.c or 02_gol_cpu_openmp_task_fort.f90)
	$(info cpu_openmp_loop : compiles both C and Fortran verisons cpu parallelised openmp codes.)
	$(info Will try compiling both loop and task based parallelism )
	$(info $(TAB)Sources : 02_gol_cpu_openmp_loop.c or 02_gol_cpu_openmp_loop_fort.f90)
	$(info )
	$(info gpu_openmp : compiles both C and Fortran verisons gpu parallelised openmp codes.)
	$(info gpu_openacc : compiles both C and Fortran verisons gpu parallelised openacc codes.)
	$(info ----------------------------------------)
	$(info NOTE:: You can append to the make argument _cc or _fort to only compile the c or fortran versions.)
	$(info For example)
	$(info cpu_serial_cc : compiles ONLY C version of serial code. )
	$(info ========================================)
	$(info )

dirs :
	[ -d obj ] || mkdir obj
	[ -d bin ] || mkdir bin

clean :
	rm obj/*
	rm bin/*

# just make an easier make name to remember
cpu_serial : cpu_serial_cc
cpu_serial_cc : bin/01_gol_cpu_serial
cpu_serial_expanded : cpu_serial_expanded_cc
cpu_serial_expanded_cc : bin/01_gol_cpu_serial_expanded
cpu_openmp : cpu_openmp_loop cpu_openmp_task
cpu_openmp_loop : cpu_openmp_loop_cc
cpu_openmp_task : cpu_openmp_task_cc
cpu_openmp_loop_cc : bin/02_gol_cpu_openmp_loop
cpu_openmp_task_cc : bin/02_gol_cpu_openmp_task

# gpu related 
gpu_openmp : bin/02_gol_gpu_openmp
gpu_openacc : bin/02_gol_gpu_openacc
gpu_mono_hip : bin/02_gol_hip_monogpu 

obj/common.o : src/common.h src/common.c
	$(CC) $(COMMONFLAGS) $(CFLAGS) -c src/common.c -o obj/common.o

bin/01_gol_cpu_serial : src/01_gol_cpu_serial.c obj/common.o
	$(CC) $(COMMONFLAGS) $(CFLAGS) -c src/01_gol_cpu_serial.c -o obj/01_gol_cpu_serial.o
	$(CC) $(COMMONFLAGS) $(CFLAGS) -o bin/01_gol_cpu_serial obj/01_gol_cpu_serial.o obj/common.o

bin/01_gol_cpu_serial_expanded : src/01_gol_cpu_serial_expanded.c obj/common.o
	$(CC) $(COMMONFLAGS) $(CFLAGS) -c src/01_gol_cpu_serial_expanded.c -o obj/01_gol_cpu_serial_expanded.o
	$(CC) $(COMMONFLAGS) $(CFLAGS) -o bin/01_gol_cpu_serial_expanded obj/01_gol_cpu_serial_expanded.o obj/common.o

bin/01_gol_cpu_opt_serial : src/01_gol_cpu_opt_serial.c obj/common.o
	$(CC) $(COMMONFLAGS) $(CFLAGS) -c src/01_gol_cpu_opt_serial.c -o obj/01_gol_cpu_opt_serial.o
	$(CC) $(COMMONFLAGS) $(CFLAGS) -o bin/01_gol_cpu_opt_serial obj/01_gol_cpu_opt_serial.o obj/common.o

bin/02_gol_cpu_openmp_loop : src/02_gol_cpu_openmp_loop.c obj/common.o
	$(OMPCC) $(OMP_FLAGS) $(COMMONFLAGS) $(CFLAGS) -c src/02_gol_cpu_openmp_loop.c -o obj/02_gol_cpu_openmp_loop.o
	$(OMPCC) $(OMP_FLAGS) $(COMMONFLAGS) $(CFLAGS) -o bin/02_gol_cpu_openmp_loop obj/02_gol_cpu_openmp_loop.o obj/common.o

bin/02_gol_cpu_openmp_task : src/02_gol_cpu_openmp_task.c obj/common.o
	$(OMPCC) $(OMP_FLAGS) $(COMMONFLAGS) $(CFLAGS) -c src/02_gol_cpu_openmp_task.c -o obj/02_gol_cpu_openmp_task.o
	$(OMPCC) $(OMP_FLAGS) $(COMMONFLAGS) $(CFLAGS) -o bin/02_gol_cpu_openmp_task obj/02_gol_cpu_openmp_task.o obj/common.o

bin/02_gol_cpu_openacc : src/02_gol_cpu_openacc.c obj/common.o
	$(OACCCC) $(OACC_FLAGS) $(COMMONFLAGS) $(CFLAGS) -c src/02_gol_cpu_openacc.c -o obj/02_gol_cpu_openacc.c
	$(OACCCC) $(OACC_FLAGS) $(COMMONFLAGS) $(CFLAGS) -o bin/02_gol_cpu_openacc obj/02_gol_cpu_openacc.o obj/common.o

bin/02_gol_hip_monogpu : src/02_gol_hip_monogpu.cpp obj/common.o
	$(HCC) $(COMMONFLAGS) $(CFLAGS) -c src/02_gol_hip_monogpu.cpp -o obj/02_gol_hip_monogpu.o
	$(HCC) $(COMMONFLAGS) $(CFLAGS) -o bin/02_gol_hip_monogpu obj/02_gol_hip_monogpu.o obj/common.o

