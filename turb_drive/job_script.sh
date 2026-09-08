#!/bin/bash
#SBATCH --job-name=athenak_turb_test
#SBATCH --partition gpupreempt
#SBATCH --constraint h200&rocky9
#SBATCH --reservation rocky9
#SBATCH --qos gpupreempt
#SBATCH --nodes 1
#SBATCH --ntasks 4            # one MPI rank per GPU
#SBATCH --ntasks-per-node 4
#SBATCH --cpus-per-task 4
#SBATCH --gpus-per-task 1     # one full GPU per rank
#SBATCH --time 01:00:00

module purge
module load modules/2.4-20250724 slurm cuda/12.5.1 openmpi/cuda-4.1.8
export LD_PRELOAD=/mnt/sw/fi/cephtweaks/lib/libcephtweaks.so
export CEPHTWEAKS_LAZYIO=1
# Minimal MPI environment (matching prior working setup)
# Uncomment if you want stack traces on abort:
# export OMPI_MCA_opal_abort_print_stack=1
# Print MPI stack on abort (new OMPI var name)
export OMPI_MCA_opal_abort_print_stack=1
# Bypass UCX/vader shared memory: force ob1 over TCP
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,tcp
# If UCX is unavoidable, prefer rc/self only (commented alternative)
# export OMPI_MCA_pml=ucx
# export UCX_TLS=rc,self
# export OMPI_MCA_btl=self,tcp
# export UCX_POSIX_USE_PROC_LINK=n
# export OMPI_MCA_btl_vader_single_copy_mechanism=none

export OMPI_MCA_opal_abort_print_stack=1
export OMPI_MCA_mpi_abort_print_stack=1
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,tcp

export athenak=/mnt/home/hdas1/athenak_nyu
export build=$athenak/build_gpu

set -e
SECONDS=0

# bash -c "unset CUDA_VISIBLE_DEVICES; \
# exec "$build/src/athena" \
# Quick sanity check: show what GPUs we actually got
# srun -N1 -n1 -l hostname nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
srun -N1 -n1 -l bash -c 'hostname; nvidia-smi --query-gpu=name,memory.total --format=csv,noheader'
 

# Use one GPU per rank and bind CPU cores accordingly
srun --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores --gpu-bind=single:1 \
  bash -c "unset CUDA_VISIBLE_DEVICES; \
  exec \"$build/src/athena\" \
  -i \"$athenak/turb_drive/turb.athinput\" \
  -d \"$athenak/turb_drive\""

echo "Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "Boom!"
