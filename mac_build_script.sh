
export athenak=/Users/hitesh/hitesh/git/athenak_nyu
export build=$athenak/bin

export PATH="$(brew --prefix llvm)/bin:$PATH"
export LDFLAGS="-L$(brew --prefix llvm)/lib $LDFLAGS"
export CPPFLAGS="-I$(brew --prefix llvm)/include $CPPFLAGS"

mkdir -p "$build"

cmake -S "$athenak" -B "$build" \
  -D CMAKE_CXX_COMPILER=clang++ \
  -D Athena_ENABLE_OPENMP=ON \
  -D Kokkos_ENABLE_OPENMP=ON \
  -D PROBLEM=turb \

cd "$build"
make -j 8