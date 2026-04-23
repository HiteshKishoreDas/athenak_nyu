#!/bin/bash

# Exit on any error
set -e

# Create log file with timestamp
log_file="build_frontier_cgm_$(date +%Y%m%d_%H%M%S).log"
echo "Build log will be saved to: ${log_file}"

# Function to echo to both stdout and log file
log_echo() {
    echo "$@" | tee -a "${log_file}"
}

# Redirect all output to both terminal and log file
exec > >(tee -a "${log_file}")
exec 2>&1

# Frontier-optimized build script for AthenaK CGM runs
log_echo "=== Building AthenaK on Frontier for CGM ==="
log_echo "=== Build started at $(date) ==="

# Load the same compiler/MPI stack expected at runtime
module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module load cmake
module unload darshan-runtime

echo "=== Loaded modules ==="
module -t list

# Set environment variables
export LD_LIBRARY_PATH="${ROCM_PATH}/lib/llvm/lib:${ROCM_PATH}/lib:${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}"

# Preserve the active runtime search path inside the executable so job scripts
# that sanitize modules can still find ROCm and Cray runtime libraries.
rpath_list="${ROCM_PATH}/lib;${ROCM_PATH}/lib/llvm/lib"
if [ -n "${CRAY_LD_LIBRARY_PATH:-}" ]; then
    rpath_list="${rpath_list};$(printf '%s' "${CRAY_LD_LIBRARY_PATH}" | tr ':' ';')"
fi

# Define paths
athenak_dir='/lustre/orion/proj-shared/ast207/hitesh/athenak_nyu'
build_dir="${athenak_dir}/build_cgm_twobin"

# Clean and create build directory
echo "=== Setting up build directory ==="

if [ -d "${build_dir}" ]; then
    echo "Removing existing build directory..."
    rm -rf "${build_dir}"
fi
mkdir -p "${build_dir}"

# Configure with CMake
echo "=== Configuring with CMake ==="
cd "${athenak_dir}"

cmake -B"${build_dir}" \
      -DAthena_ENABLE_MPI=ON \
      -DKokkos_ARCH_ZEN3=ON \
      -DKokkos_ARCH_VEGA90A=ON \
      -DKokkos_ENABLE_HIP=ON \
      -DCMAKE_CXX_COMPILER=CC \
      -DCMAKE_CXX_FLAGS="-I${ROCM_PATH}/include -munsafe-fp-atomics" \
      -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/lib -L${ROCM_PATH}/lib/llvm/lib -lamdhip64 -Wl,-rpath,${ROCM_PATH}/lib -Wl,-rpath,${ROCM_PATH}/lib/llvm/lib" \
      -DCMAKE_BUILD_RPATH="${rpath_list}" \
      -DCMAKE_INSTALL_RPATH="${rpath_list}" \
      -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON \
      -DPROBLEM=cgm_cooling_flow_full

# Build
echo "=== Building AthenaK ==="
cd "${build_dir}"
make -j16

echo "=== Build complete! ==="
echo "Executable location: ${build_dir}/src/athena"
echo "=== Build finished at $(date) ==="
echo ""
echo "Build log saved to: ${log_file}"
