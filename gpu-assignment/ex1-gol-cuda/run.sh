#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --account=courses0100
#SBATCH --nodes=1
#SBATCH --partition=gpuq

module use /group/courses0100/software/nvhpc/modulefiles
module load nvhpc/21.9
module load gcc/9.2

DIR="$( pwd )"
echo "The current script directory is: $DIR"
echo 

echo "Make Directories"
[ -d obj ] || mkdir obj
[ -d bin ] || mkdir bin
echo


echo "Compiling code"
srun --export=all -u -n 1 make clean
srun --export=all -u -n 1 make
echo

# Run C code
# echo "Running C code"
# srun --export=all -u -n 1 bin/game_of_life_c 10 10 5
# echo 

echo "Running C CUDA"
srun --export=all -u -n 1 bin/game_of_life_cuda_c 1000 10000 10
# srun --export=all u -n 1 make -f bin/game_of_life_c
# srun --export=all u -n 1 make -f bin/game_of_life_cuda_c
``