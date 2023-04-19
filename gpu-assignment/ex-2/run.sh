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

Run C code
echo "Running C code"
srun --export=all -u -n 1 bin/01_gol_cpu_serial 10000 10000 100
echo

# Run Directive code
echo "Running Directive code"
srun --export=all -u -n 1 bin/02_gol_gpu_openacc 10000 10000 100
echo

echo "Generating profiling report for C"
srun --export=all -u -n 1 gprof -lbp bin/02_gol_gpu_openacc gmon.out
echo


