# Non-Equilibrium Chemical Network Test Suite

## Overview

`test_nonequilibrium.py` provides a comprehensive testing framework for validating the chemical network ODE solver, Jacobian computations, and physical constraint enforcement.

## Running Tests

```bash
# Run all tests
python src/chemistry/test_nonequilibrium.py

# Run with verbose output
python src/chemistry/test_nonequilibrium.py  # Shows rates and Jacobian spectrum by default
```

## Test Categories

### Unit Tests

#### Test 1: Rate Function Evaluation ✅
- **What it does**: Verifies that the rate function correctly evaluates reaction rates and RHS derivatives
- **Checks**:
  - Output shapes match expected dimensions
  - No NaN or Inf values
  - Most rates are physically reasonable (positive)
- **Output**: `[PASSED] Max rate, Min rate, Positive rates count`

#### Test 2: Stoichiometric Matrix ✅
- **What it does**: Validates the stoichiometric matrix construction
- **Checks**:
  - Shape: (n_species, n_reactions)
  - All elements are integers (valid stoichiometry)
  - Matrix structure is consistent
- **Output**: `[PASSED] Shape, Species/Reaction counts`

#### Test 3: Analytic Jacobian Computation ✅
- **What it does**: Verifies analytic Jacobian computation from reaction stoichiometry
- **Checks**:
  - Output shape: (n_species, n_species)
  - No NaN or Inf values
  - Spectral properties (eigenvalues, condition number)
  - Stiffness ratio (important for ODE solver behavior)
- **Output**:
  ```
  [PASSED] 
    Max eigenvalue (real): X.XXe+YY
    Min eigenvalue (real): -X.XXe+YY
    Condition number: X.XXe+YY
    Stiffness ratio: X.XXe+YY
  ```

### Integration Tests

#### Test 4: Atom Conservation Enforcement ✅
- **What it does**: Tests the atom conservation projection algorithm
- **Checks**:
  - Conservation enforced for all elements (H, He, C, N, O, etc.)
  - Relative conservation error < 1e-6
- **Output**: `[PASSED] Element count, Conservation errors count`

#### Test 5: SUNDIALS (CVODE) Integration ✅
- **What it does**: Runs full ODE integration using SUNDIALS CVODE solver
- **Checks**:
  - Integration completes successfully
  - Final state is physically reasonable
  - Analytic Jacobian is accepted by solver
  - Output: `[CVODE] analytic Jacobian callback passed to solver`
- **Output**:
  ```
  [PASSED]
    Max concentration: X.XXe-YY
    Min concentration: X.XXe-YY
    Non-zero species: N/20
  ```

### Non-Equilibrium Dynamics Tests

#### Test 6: Non-Equilibrium Evolution ✅
- **What it does**: Verifies temperature-dependent chemistry evolution
- **Temperature range**: 10K, 30K, 100K, 300K, 1000K
- **Checks**:
  - Chemistry changes significantly with temperature
  - Rate variation coefficient > 0.1 (indicates stiff system)
  - Active reactions vary with temperature
- **Output**:
  ```
  [PASSED]
    Temperature range: 10.0K - 1000.0K
    Rate variation coefficient: X.XXX (significant: True)
  ```

## Configuration

Modify `TestConfig` in the test file to customize:

```python
config = TestConfig()
config.temperatures = [10.0, 30.0, 100.0, 300.0, 1000.0]  # Temperature range
config.t_end = 1e3              # Integration time (seconds)
config.rtol = 1e-6              # Relative tolerance
config.atol = 1e-12             # Absolute tolerance
config.verbose = True           # Verbose output
config.show_rates = True        # Show reaction rate statistics
config.show_jacobian_spectrum = True  # Show eigenvalue analysis
```

## Output Interpretation

### Success Indicators
- ✅ All 6 tests pass
- `[CVODE] analytic Jacobian callback passed to solver` - Analytic derivatives are used
- Condition number ~1e30 - System is appropriately stiff (typical for chemistry)
- Rate variation coefficient > 0.1 - Non-equilibrium dynamics present

### Warning Signs
- ❌ Test failure - Check error message for specific issue
- High NaN/Inf count - Integration may be diverging
- Perfect atom conservation - May indicate trivial test case
- Low rate variation - Chemistry may not be temperature-dependent

## Physical Interpretation

### Stiffness Ratio
- Ratio = |max eigenvalue| / |min eigenvalue|
- High ratio (> 1e10) → Stiff system (requires implicit solver)
- Low ratio (< 10) → Non-stiff system (can use explicit solver)
- Typical chemistry: Very high ratio (1e15+) → CVODE essential

### Eigenvalue Spectrum
- **Positive eigenvalues** → Unstable modes (energy injection)
- **Negative eigenvalues** → Stable modes (energy dissipation)
- **Complex eigenvalues** → Oscillatory behavior
- **Small magnitude** → Nearly conserved species (slow dynamics)

### Atom Conservation Error
- Should be < 1e-10 after projection
- Indicates numerical integration accuracy
- Larger errors suggest need for tighter tolerances

## Extending Tests

To add custom checks, extend the functional tests (e.g., in `test_nonequilibrium.py`) with new helper functions that run your assertions and append to the `results` dict.

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `SUNDIALS not available` | Install scikits.odes: `pip install scikits.odes` |
| Test 5 skipped | SUNDIALS not installed; SciPy will be tested instead |
| High NaN count | Reduce integration time or tighten tolerances |
| Atom conservation failed | Check species elemental composition data |
| Jacobian condition number very high | This is expected for stiff chemistry systems |

## Dependencies

- `numpy` (optional, but recommended for performance)
- `scipy` (required for solve_ivp fallback)
- `scikits.odes` (optional, for SUNDIALS/CVODE)

## Exit Codes

- `0` - All tests passed
- `1` - One or more tests failed

## Performance Notes

- Full test suite: ~10-30 seconds (depends on CVODE availability)
- Each CVODE integration: ~1-5 seconds
- Jacobian computation: ~milliseconds per evaluation
- Test suite memory: ~10-50 MB

---

For full system documentation, see `DOCUMENTATION.md`
