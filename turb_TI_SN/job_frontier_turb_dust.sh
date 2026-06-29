#!/bin/bash
#SBATCH -A AST207
#SBATCH -J dust_TI
#SBATCH -o logs/%x.%j.out
#SBATCH -t 02:00:00
#SBATCH -p batch
#SBATCH --qos=debug
#SBATCH -N 1

# Load modules
module restore
module load cpe/24.07 PrgEnv-amd cray-mpich/8.1.30 craype-accel-amd-gfx90a amd/6.2.0 rocm/6.2.0
module load cmake cray-python emacs
module unload darshan-runtime
module -t list

# Environment setup
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export FI_MR_CACHE_MONITOR=disable
export MPICH_SMP_SINGLE_COPY_MODE=NONE

# Paths
run_dir=/lustre/orion/proj-shared/ast207/hitesh/athenak_nyu/build_turb/src/
executable=athena
input_dir=/lustre/orion/proj-shared/ast207/hitesh/athenak_nyu/turb_TI_SN/
output_dir=/lustre/orion/proj-shared/ast207/hitesh/athenak_nyu/turb_TI_SN/
input_file=turb_dust.athinput

cd $run_dir
mkdir -p ${output_dir}
mkdir -p ${output_dir}/bin
mkdir -p ${output_dir}/rst
mkdir -p logs

# Check for restart files
test_file=$(find $output_dir/rst -maxdepth 1 -name "*.rst" -print -quit)
if [ -n "$test_file" ]; then
  restart_files=$(ls -t $output_dir/rst/*.rst)
  restart_file=(${restart_files[0]})
  restart_line="-r $restart_file -i ${input_dir}/${input_file}"
  printf "\nrestarting $name from $restart_file\n\n"
else
  restart_line="-i ${input_dir}/${input_file}"
  printf "\nstarting $name from beginning\n\n"
fi

# Run the simulation
srun -N 1 -n 8 --ntasks-per-node=8 \
     --gpus-per-task=1 --gpu-bind=closest \
     ./${executable} \
     $restart_line \
     -t 02:20:00 \
     -d ${output_dir} \
     $arguments

echo "Job completed at $(date)"

# Print performance information
printf "\n========== Performance Summary ==========\n"
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Elapsed,State
printf "=========================================\n"
