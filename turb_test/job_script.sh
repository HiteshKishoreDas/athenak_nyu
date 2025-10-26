#!/bin/bash
#SBATCH --job-name=athenak_turb_test
#SBATCH --partition gpu
#SBATCH --constraint a100-80gb
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 16
#SBATCH --gpus-per-task 1
#SBATCH --time 01:00:00

module purge
module load modules/2.4-20250724 slurm cuda/12.5.1 openmpi/cuda-4.1.8
export LD_PRELOAD=/mnt/sw/fi/cephtweaks/lib/libcephtweaks.so
export CEPHTWEAKS_LAZYIO=1

export athenak=/mnt/home/hdas1/athenak_nyu
export build=$athenak/build

srun --cpus-per-task=$SLURM_CPUS_PER_TASK --cpu-bind=cores --gpu-bind=single:2 \
bash -c "unset CUDA_VISIBLE_DEVICES"; \
$build/src/athena -i $athenak/turb_test/convection.athinput -d $athenak/turb_test
