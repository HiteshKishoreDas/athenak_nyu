# Quick Start: Running the Non-Equilibrium Test Suite

## One-Liner to Run Tests

```bash
cd /Users/hitesh/hitesh/git/athenak_nyu && python src/chemistry/test_nonequilibrium.py
```

## Expected Success Output

```
==================================================================
====  CHEMICAL NETWORK NON-EQUILIBRIUM TEST SUITE
==================================================================
====
✅ TEST 1 (Rate Function): PASSED at T=100.0K
   Max rate: X.XXe-YY, Min rate: X.XXe-YY
   Positive rates: 67/67

✅ TEST 2 (Stoichiometry): PASSED
   Shape: (20 species, 67 reactions)

✅ TEST 3 (Jacobian): PASSED at T=100.0K
   Shape: (20, 20)
   Condition number: 5.082e+34
   Stiffness ratio: 1.237e-17

✅ TEST 4 (Atom Conservation): PASSED
   Elements: 22, Conservation errors: 0

✅ TEST 5 (SUNDIALS Integration): PASSED at T=100.0K
   [CVODE] analytic Jacobian callback passed to solver
   Max concentration: 1.700e-07
   Non-zero species: 18/20

✅ TEST 6 (Non-equilibrium Evolution): PASSED
   Temperature range: 10.0K - 1000.0K

==================================================================
==== TEST SUMMARY
==================================================================
  test_rate_function                       ✅ PASSED
  test_stoichiometric_matrix               ✅ PASSED
  test_jacobian_computation                ✅ PASSED
  test_atom_conservation                   ✅ PASSED
  test_integration_sundials                ✅ PASSED
  test_nonequilibrium_evolution            ✅ PASSED

Total: 6/6 tests passed
==================================================================
✅ All tests passed!
```

## What Each Test Does

| # | Test | Time | What It Validates |
|---|------|------|-------------------|
| 1 | Rate Function | <1s | Chemical reaction rates computed correctly |
| 2 | Stoichiometry | <1s | Stoichiometric matrix structure is valid |
| 3 | Jacobian | <1s | Exact partial derivatives computed correctly |
| 4 | Atom Conservation | <1s | Element conservation enforcement works |
| 5 | SUNDIALS Integration | ~20s | ODE solver with analytic Jacobian works |
| 6 | Non-equilibrium | ~7s | Temperature-dependent chemistry captured |

## Test Meanings

**✅ PASSED** = All assertions passed, system working correctly

**⏭️ SKIPPED** = Test skipped (usually SUNDIALS not available, will use SciPy instead)

**❌ FAILED** = Test failed, see error message for details

## For Custom Configurations

Edit the `TestConfig()` parameters:

```python
# In src/chemistry/test_nonequilibrium.py, line ~540
config = TestConfig()
config.temperatures = [10.0, 50.0, 100.0, 500.0, 1000.0]  # Your temps
config.t_end = 5e3  # Integration time (seconds)
config.verbose = True
config.show_jacobian_spectrum = True  # See eigenvalue analysis
test_suite = ChemicalNetworkTest(network, config)
results = test_suite.run_all_tests()
```

## Typical Runtime

- **Without SUNDIALS**: ~5 seconds (uses SciPy)
- **With SUNDIALS**: ~30 seconds (full CVODE integration)

## Key Files

| File | Purpose |
|------|---------|
| `test_nonequilibrium.py` | Main test suite |
| `TEST_README.md` | Detailed test documentation |
| `DOCUMENTATION.md` | System architecture |
| `TEST_SETUP_SUMMARY.md` | Setup and extension guide |

## Dependencies

```bash
# Required
pip install numpy scipy

# Optional but recommended
pip install scikits.odes  # For SUNDIALS/CVODE
```

## Exit Codes

- `0` = All tests passed ✅
- `1` = One or more tests failed ❌

---

**That's it!** The test suite is ready to validate your chemical network solver.
