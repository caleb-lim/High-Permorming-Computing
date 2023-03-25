#!/bin/bash --login 
#SBATCH --tasks=1 
#SBATCH --gres=gpu:1
#SBATCH --mem=1G
#SBATCH --partition=gpuq 
#SBATCH --time=00:05:00 
#SBATCH --export=NONE
module purge
module load cascadelake
module load pgi/19.7
module load gcc/10.2
module load cuda

# Compile
echo "Compiling code"
srun --export=all -u -n 1 make clean
srun --export=all -u -n 1 make
echo

# Run C code
echo "Running C code"
srun --export=all -u -n 1 bin/01_gol_cpu_serial 10 10 5
echo

# Run Directive code
echo "Running Directive code"
srun --export=all -u -n 1 bin/02_gol_gpu_openacc 10 10 5
