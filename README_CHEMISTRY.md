Chemistry ODE demo

This repository contains a small chemistry network loader and ODE-system builder.

Quick setup (macOS, using the project's virtualenv):

1) Ensure Homebrew packages (only if you will build scikits.odes from source):

   brew install sundials open-mpi

2) Activate your venv and install Python deps:

   source .venv/bin/activate
   pip install -r requirements.txt

   If scikits.odes fails to build, install SUNDIALS and Open MPI (see step 1),
   then re-run the pip install with environment variables:

   CFLAGS="-I$(brew --prefix sundials)/include -I$(brew --prefix)/include" \
   LDFLAGS="-L$(brew --prefix sundials)/lib -L$(brew --prefix)/lib" \
   PKG_CONFIG_PATH="$(brew --prefix sundials)/lib/pkgconfig" \
   pip install scikits.odes

3) Run the ODE demo (uses CVODE if available, otherwise falls back to SciPy):

   .venv/bin/python src/chemistry/build_ode_system.py

Files of interest:
- `src/chemistry/custom_network_loader.py` : KIDA file parser
- `src/chemistry/build_ode_system.py` : builds stoichiometric matrix and rate function, demo integration
- `src/chemistry/custom/` : KIDA custom files used in the demo
