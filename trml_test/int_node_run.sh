#!/bin/bash

module purge
module load modules/2.4-20250724 slurm cuda/12.5.1 openmpi/cuda-4.1.8
export LD_PRELOAD=/mnt/sw/fi/cephtweaks/lib/libcephtweaks.so
export CEPHTWEAKS_LAZYIO=1

export athenak=/mnt/home/hdas1/athenak_nyu
export build=$athenak/build

set -e
SECONDS=0

srun --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores --gpu-bind=single:2 \
  bash -c "unset CUDA_VISIBLE_DEVICES; \
  exec \"$build/src/athena\" \
  -i \"$athenak/trml_test/convection.athinput\" \
  -d \"$athenak/turb_test\""



echo "Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "Boom!"
