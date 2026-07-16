module purge
module load modules/2.4-20250724 slurm cuda/12.5.1 openmpi/cuda-4.1.8
# export LD_PRELOAD=/mnt/sw/fi/cephtweaks/lib/libcephtweaks.so
# export CEPHTWEAKS_LAZYIO=1
# module load modules/2.4-20250724 slurm cuda/11.8.0 openmpi/cuda-4.0.7
export LD_PRELOAD=/mnt/sw/fi/cephtweaks/lib/libcephtweaks.so
export CEPHTWEAKS_LAZYIO=1

# build=/mnt/home/hdas1/athenak_nyu/bin
# athenak=/mnt/home/hdas1/athenak_nyu

athenak=/mnt/home/hdas1/athenak_nyu
build=$athenak/build_gpu

mkdir -p "$build"

cmake -S "$athenak" -B "$build" \
  -D CMAKE_CXX_COMPILER="$athenak/kokkos/bin/nvcc_wrapper" \
  -D Kokkos_ENABLE_CUDA=On \
  -D Kokkos_ARCH_AMPERE80=On \
  -D Athena_ENABLE_MPI=On \
  -D PROBLEM=turb
  # -D CMAKE_BUILD_TYPE=Debug \
  # -D Kokkos_ENABLE_DEBUG=On

cd $build
make -j 16