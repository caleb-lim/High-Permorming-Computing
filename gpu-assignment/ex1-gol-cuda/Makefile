.PHONY : dirs

COMMONFLAGS=-O0 -g#-DUSEPNG
INCLUDES=#-I/usr/include
LIBS=#-L/usr/lib64/ -lpng
NVFORTRAN_FLAGS=-module obj

all : dirs bin/game_of_life_c bin/game_of_life_cuda_c

dirs :
	[ -d obj ] || mkdir obj
	[ -d bin ] || mkdir bin

clean :
	rm obj/*
	rm bin/*

# ============== BEGIN CUDA C++ =======================================

obj/common.o : src/common.h src/common.c
	$(CC) $(COMMONFLAGS) $(INCLUDES) -o $@ -c src/common.c $(LIBS)

bin/game_of_life_c : src/game_of_life.c obj/common.o
	gcc $(COMMONFLAGS) $(INCLUDES) -o $@ src/game_of_life.c obj/common.o $(LIBS)

bin/game_of_life_cuda_c : src/game_of_life.cu obj/common.o
	nvcc -G -lineinfo $(COMMONFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

# ============== END CUDA C++ ===========================================

# ============== BEGIN CUDA FORTRAN =====================================

# obj/cpu_game_of_life_mod.o: src/game_of_life_mod.f90
# 	nvfortran $(GFORTRAN_FLAGS) -o $@ -c $<

# bin/game_of_life_f : src/game_of_life.f90 obj/cpu_game_of_life_mod.o
# 	nvfortran $(GFORTRAN_FLAGS) $(FLAGS) -o $@ $^

# obj/gpu_game_of_life_mod.o : src/game_of_life_mod.cuf
# 	nvfortran $(NVFORTRAN_FLAGS) -o $@ -c $<

# bin/game_of_life_cuda_f : src/game_of_life.cuf obj/gpu_game_of_life_mod.o obj/cpu_game_of_life_mod.o 
# 	nvfortran $(NVFORTRAN_FLAGS) $(COMMONFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

