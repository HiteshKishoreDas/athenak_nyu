# Chemical Network ODE Solver Documentation

## Overview

This chemistry module provides a complete pipeline for loading KIDA (Kinetic Database for Astrochemistry) chemical networks and solving the resulting ordinary differential equations (ODEs) using both SUNDIALS (CVODE) and SciPy integrators with analytic Jacobian support.

### Key Features

- **Custom KIDA Network Loading**: Parse species, reactions, surface data, and metadata from KIDA files
- **Stoichiometric Matrix Construction**: Build the matrix relating species changes to reaction rates
- **Analytic Jacobian Computation**: Compute exact derivatives from reaction stoichiometry for improved solver efficiency
- **Dual Integrator Support**: SUNDIALS (CVODE) for production code, SciPy for fallback/validation
- **Physical Constraints**: Enforce atom conservation and non-negativity after integration
- **ANSI-aware Output**: Pretty-printed tables with color coding for problematic values

---

## Module Architecture

```
custom_network_loader.py  ──→  load_network() (loads KIDA files into a dict)
                                    ↓
build_ode_system.py       ──→  make_ode_system()
                               (constructs ODE system)
                                    ↓
                          ┌──────────┼──────────┐
                          ↓          ↓          ↓
                    rate_func    species     S (stoichiometric)
                    (RHS eval)    (list)      matrix
                          ↓          ↓          ↓
                    integrate_with_sundials()  /
                    run_demo() ──────────────────/
                          ↓
                    jacobian_utils.py
                    (analytic_jacobian)
                          ↓
                    cli_utils.py
                    (print_clamped_table,
                     enforce_atom_conservation)
```

---

## Detailed Step-by-Step Process

### Step 1: Load Chemical Network
**File: `custom_network_loader.py`**

#### What happens:
`load_network(directory)` reads all KIDA files in a directory and returns a plain dict with `species`, `reactions`, `element_order`, `metadata`, and `surface` entries.

#### Input files:
- `kida_spec_*.dat` - Species definitions (name, charge, element composition)
- `kida_spec_readme_*.txt` - Species metadata and element order
- `kida_reac_*.dat` - Gas-phase reactions
- `kida_reac_tb_*.dat` - Additional reaction sets
- `kida_surface_*.dat` - Surface processes (barriers, diffusion, desorption, branching ratios)
- `kida_readme_*.txt` - Metadata about rate coefficients

#### Key data structures created (dict form):
```python
network["species"] = {
    "H": {"charge": 0, "composition": {"H": 1}},
    "H+": {"charge": 1, "composition": {"H": 1}},
    "H2": {"charge": 0, "composition": {"H": 2}},
    # ... more species
}

network["reactions"] = [
    {
        "reactants": ["H", "H", "H"],
        "products": ["H2", "H"],
        "alpha": 0.0,
        "beta": -0.5,
        "gamma": 1e-31,
        "Tmin": 10.0,
        "Tmax": 10000.0,
        "types": [0],
    },
    # ... more reactions
]

network.element_order = ["H", "He", "C", "N", "O", ...]  # Determines atom conservation order
```

#### Key methods:
- `load_from_directory(path)` - Main entry point
- `_load_species(file)` - Parse species file (format: name charge H He C N O ...)
- `_load_reactions(file)` - Parse reaction file
- `calculate_rate_coefficient(rxn, T, av, cr_zeta)` - Compute k(T) using Arrhenius/power-law formula

---

### Step 2: Build ODE System
**File: `build_ode_system.py` → `make_ode_system()`**

#### What happens:
Constructs the mathematical ODE system: the stoichiometric matrix **S** and the rate evaluation function.

#### Stoichiometric Matrix Construction:
For each species and reaction, determine how many molecules are consumed/produced.

**Example:**
```
Reaction: H + H + H → H2 + H

Species order: [H, H+, H-, H2, H2+, ...]
              [0,  1,   2,  3,   4,  ...]

For reaction j, stoichiometric column S[:, j]:
  S[0, j] = -2   (2 H consumed)
  S[3, j] = +1   (1 H2 produced)
  S[1, j] = 0    (H+ unaffected)
  ... (rest are 0)
```

#### Reaction Data Storage:
For each reaction, store:
```python
reaction_data[j] = {
    "reactants": ["H", "H", "H"],
    "products": ["H2", "H"],
    "reac_inds": [0, 0, 0],          # Species indices in `species` list
    "reac_mult": [1, 1, 1],          # Multiplicities (powers in rate law)
    "reac_counts": {"H": 3},         # Element multiplicity counts
    "prod_counts": {"H2": 1, "H": 1},
    "rxn": {...},                    # Original reaction dict
}
```

#### Rate Function Signature:
```python
rates, rhs = rate_func(concs, t, T, av, cr_zeta)
```

**Inputs:**
- `concs` - 1D array of species concentrations (length = n_species)
- `t` - Current time (for ODE signature compatibility)
- `T` - Temperature (K)
- `av` - Visual extinction (for photorates)
- `cr_zeta` - Cosmic-ray ionization rate (s⁻¹)

**Outputs:**
- `rates` - 1D array of reaction rates (length = n_reactions)
- `rhs` - 1D array of d[species]/dt (length = n_species), computed as: **rhs = S · rates**

#### Mass-Action Kinetics:
For reaction k with rate coefficient k(T) and reactants with concentrations c and multiplicities n:
```
rate_k = k(T) × ∏(c_j^n_j)  for all reactants j
```

#### Returns:
```python
species, S, rate_func, reaction_data = make_ode_system(network)
```

---

### Step 3: Compute Analytic Jacobian
**File: `jacobian_utils.py` → `analytic_jacobian()`**

#### What happens:
Computes the exact Jacobian matrix J ∈ ℝ^(n_species × n_species) from reaction stoichiometry.

#### Mathematical Foundation:
For a system of ODEs: d[species_i]/dt = (RHS)_i

The Jacobian is: J[i,j] = ∂(RHS)_i / ∂c_j

Using the stoichiometric relationship:
```
(RHS)_i = Σ_k S[i,k] × rate_k

Therefore:
J[i,j] = Σ_k S[i,k] × (∂rate_k / ∂c_j)
```

#### Partial Derivative Computation:
For mass-action rate: rate_k = k(T) × ∏(c_l^n_l) for reactant l

The partial derivative with respect to reactant j with multiplicity n_j:
```
∂rate_k / ∂c_j = (rate_k / c_j) × n_j   (if j is a reactant)
                = 0                        (if j is NOT a reactant)
```

Handle c_j → 0 by using: ∂rate_k / ∂c_j ≈ k(T) × n_j × c_j^(n_j - 1)

#### Algorithm:
```python
1. Evaluate all reaction rates: rate_k = k(T) × ∏(c_l^n_l)
2. Initialize J = 0
3. For each reaction k:
4.   For each reactant j with multiplicity n_j:
5.     Compute dr_k/dc_j = (rate_k / c_j) × n_j
6.     For each species i:
7.       J[i,j] += S[i,k] × dr_k/dc_j
8. Return J (dense matrix)
```

#### Implementation Details:
- **Numpy version** (preferred): Uses vectorization for speed
- **Pure Python fallback** (if numpy unavailable): Uses loops
- Handles small concentrations (c_j ≈ 0) by checking for division by zero

#### Why Analytic Jacobian?
- **Accuracy**: Exact derivatives instead of finite-difference approximation
- **Efficiency**: No extra ODE function evaluations needed
- **Stiffness**: Improves solver robustness for stiff systems (common in chemistry)

---

### Step 4: Integrate ODE System
**File: `build_ode_system.py` → `integrate_with_sundials()` and `run_demo()`**

#### Two Integrator Paths:

##### Path A: SUNDIALS (CVODE)
**When:** Available (scikits.odes installed)
**Method:** BDF (Backward Differentiation Formula) for stiff systems
**Jacobian:** Analytic (if provided), else numeric finite-difference
**Signature:**
```python
sol = integrate_with_sundials(
    rate_func, y0, t0, tf,
    T=100.0, av=10.0, cr_zeta=1.3e-17,
    reaction_data=reaction_data,
    species=species,
    S=S,
    network=network,
    use_analytic_jac=True
)
```

**Steps:**
1. Build wrapper functions for rate (RHS) and Jacobian callbacks
2. Create CVODE solver instance with callbacks
3. Integrate from t0 to tf, evaluating solution at 100 time points
4. Return result object with `.y` (concentrations at each time point)

**Output:** `[CVODE] analytic Jacobian callback passed to solver`

##### Path B: SciPy (solve_ivp)
**When:** SUNDIALS fails or unavailable
**Method:** BDF (implicit multistep)
**Jacobian:** Analytic (if provided), else estimated by solver
**Signature:**
```python
sol = solve_ivp(
    rhs_scipy, (t0, tf), x0,
    method="BDF",
    rtol=1e-6, atol=1e-12,
    jac=jac_scipy
)
```

**Steps:**
1. Wrap rate function for SciPy signature: `rhs_scipy(t, y) → dy/dt`
2. Wrap Jacobian for SciPy signature: `jac_scipy(t, y) → J`
3. Call solve_ivp with analytic Jacobian
4. Extract final concentrations from `.y[:, -1]`

#### Integration Parameters:
- `t0, tf` - Start and final times
- `T` - Temperature (constant during integration)
- `av` - Visual extinction (constant)
- `cr_zeta` - Cosmic-ray ionization rate (constant)
- `rtol=1e-6, atol=1e-12` - Relative and absolute tolerances

#### Output:
Solution object with:
- `.y` - 2D array (n_species × n_timepoints) of concentrations
- `.t` - 1D array of time points
- `.status` - Integration status (0 = success)

---

### Step 5: Enforce Physical Constraints
**File: `cli_utils.py` → `enforce_atom_conservation()` and clamping**

#### Why?
Numerical integration can produce:
- **Negative concentrations** (physically impossible)
- **Atom conservation violations** (due to roundoff errors)

#### Step 5A: Atom Conservation
**Function:** `enforce_atom_conservation(final, initial, species, network)`

**Algorithm:**
1. For each element E:
   - Compute total atoms of E at end: `n_E_final = Σ_i (composition[i][E] × conc[i])`
   - Compute total atoms of E at start: `n_E_initial = Σ_i (composition[i][E] × x0[i])`
2. Compute deficit for each element: `deficit[E] = n_E_initial - n_E_final`
3. Distribute deficit proportionally among species containing E:
   - For each species i containing E:
     - `adjustment[i] += deficit[E] / num_species_with_E`
   - `final_corrected[i] = final[i] + adjustment[i]`

**Example:**
```
Initial: H = 0.28 (all in H atoms)
Final:   H = 0.27 (lost to H2, but roundoff error)
Deficit: 0.01

Distribute 0.01 among H, H+, H-, H2, ... proportionally
```

#### Step 5B: Non-Negativity Clamping
**Algorithm:**
```python
final_clamped[i] = max(final_conserved[i], 0.0)
```

**Tracking:**
- If `final_clamped[i] != final_conserved[i]`, the value was clamped
- Store original value for reporting: `(was X.XXe-YY)`

---

### Step 6: Display Results
**File: `cli_utils.py` → `print_clamped_table()`**

#### ANSI-aware Table Printing
**Features:**
- Strips ANSI codes before computing column width (prevents misalignment)
- Colors clamped values in red: `\033[91m`
- Shows original values for clamped entries: `(was -1.23e-45)`
- Aligns columns with proper padding

**Example Output:**
```
+------------------------+--------------------------------------+
| Species                |      Clamped (cm^-3) (was ...)       |
+------------------------+--------------------------------------+
| H                      |           0.000e+00 (was -1.720e-16) |
| H+                     |                            7.000e-08 |
| H-                     |                            9.946e-24 |
| H2                     |                            4.239e-17 |
+------------------------+--------------------------------------+
```

**Red coloring:** Applied to clamped values so users can quickly spot non-physical corrections.

---

## Data Flow Diagram

```
1. KIDA Files (*.dat, *.txt)
   ↓
2. CustomChemicalNetwork.load_from_directory()
   → species dict, reactions list, metadata
   ↓
3. make_ode_system(network)
   → species[], S (matrix), rate_func, reaction_data[]
   ↓
4. prepare_initial_concentrations(species)
   → x0[] (e.g., H=2.8e-7, He=4e-8, others=1e-8)
   ↓
5. rate_func(x0, t=0, T, av, cr_zeta)
   → rates[], rhs[] (sample evaluation)
   ↓
6a. integrate_with_sundials(rate_func, x0, t0, tf, ...)
    + analytic_jacobian(concs, reaction_data, ...)
    ↓ (or fallback to SciPy)
6b. solve_ivp(rhs_scipy, (t0, tf), x0, jac=jac_scipy)
    ↓
7. final_concs[] from integrator
   ↓
8. enforce_atom_conservation(final_concs[], x0[], species, network)
   → final_conserved[]
   ↓
9. Clamp to non-negative: final_clamped[] = max(final_conserved[], 0)
   ↓
10. print_clamped_table(species, final_conserved, final_clamped)
    ↓
    Output: ANSI table with red-colored clamped values
```

---

## Key Algorithms

### A. Rate Coefficient Calculation (Temperature Dependent)
```python
# From network.calculate_rate_coefficient():
k(T) = α × (T/300)^β × exp(γ/T)  [Basic Arrhenius]
     or specific formulas for cosmic-ray, photodissociation, etc.
```

### B. Stoichiometric Matrix Product (ODE RHS)
```python
# Mass-action kinetics
rates[j] = k_j(T) × ∏(concs[i]^n_ij)  for reactants i in reaction j
RHS[i] = Σ_j S[i,j] × rates[j]        (matrix-vector product)
```

### C. Jacobian Assembly
```python
J[i,j] = ∂(RHS)_i / ∂(concs)_j
       = Σ_k S[i,k] × ∂rate_k / ∂concs_j
       = Σ_k S[i,k] × (rate_k / concs_j) × n_kj  (if j is reactant in k)
```

---

## Usage Example

```python
from chemistry.custom_network_loader import load_network
from chemistry.build_ode_system import make_ode_system, run_demo, prepare_initial_concentrations

# Step 1: Load network (dict-based)
network = load_network("src/chemistry/custom")
if not network.get("species") or not network.get("reactions"):
    raise SystemExit("Failed to load network")

# Step 2: Build ODE system
species, S, rate_func, reaction_data = make_ode_system(network)

# Step 3: Prepare initial concentrations
x0 = prepare_initial_concentrations(species)

# Step 4: Run demo (sample evaluation + integration)
run_demo(
    rate_func, x0, species, network,
    t0=0.0, tf=1e3, T=100.0, av=10.0, count=10,
    S=S, reaction_data=reaction_data
)
```

---

## Files Reference

| File | Purpose | Key Classes/Functions |
|------|---------|----------------------|
| `custom_network_loader.py` | Load KIDA files | `load_network()`, `calculate_rate_coefficient()` |
| `build_ode_system.py` | ODE system construction & integration | `make_ode_system()`, `integrate_with_sundials()`, `run_demo()` |
| `jacobian_utils.py` | Analytic Jacobian computation | `analytic_jacobian()` |
| `cli_utils.py` | Output formatting & constraints | `print_clamped_table()`, `enforce_atom_conservation()` |

---

## Performance Considerations

1. **Jacobian Selection:**
   - Analytic > Numeric FD > None
   - Analytic Jacobian saves ~(n_species + 1) ODE function evaluations per Newton iteration
   - Critical for stiff systems with many species

2. **Tolerance Settings:**
   - `rtol=1e-6, atol=1e-12` - standard chemistry
   - Lower (stricter) for smaller concentrations
   - Tighter tolerances → more Newton iterations → longer runtime

3. **Time Integration:**
   - BDF (implicit) preferred over RK45 (explicit) for chemical networks
   - Implicit methods are A-stable, handle stiffness well
   - SUNDIALS (CVODE) faster than SciPy for production use

---

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| `NameError: name 'reaction_data' is not defined` | Missing parameter in function call | Ensure `S` and `reaction_data` passed to `run_demo()` |
| Negative concentrations in output | Numerical roundoff | Clamping applied; check if issue is physical |
| Atom conservation violated | Stiff solver error | Tighten `rtol` / `atol` or use analytic Jacobian |
| SUNDIALS integration fails | Missing scikits.odes | Falls back to SciPy; install SUNDIALS bindings if needed |
| Jacobian callback not used | Solver doesn't accept it | Try SciPy fallback; check CVODE build flags |

---

## Future Extensions

- [ ] Implicit surface chemistry (Langmuir-Hinshelwood kinetics)
- [ ] Temperature evolution (coupled to thermal ODE)
- [ ] Sparse Jacobian support (for large networks >100 species)
- [ ] Adaptive time-stepping with error control
- [ ] Parallel reaction rate evaluation (GPU/OpenMP)
- [ ] Reaction network visualization
- [ ] Sensitivity analysis (∂RHS / ∂rate_coefficient)
