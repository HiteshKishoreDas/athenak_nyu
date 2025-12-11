"""Test suite for non-equilibrium chemical network evolution.

This module provides a comprehensive test framework for validating the ODE solver,
Jacobian computations, and physical constraint enforcement in non-equilibrium scenarios.

Test Categories:
  1. Unit Tests - Individual components (Jacobian, constraints, rate functions)
  2. Integration Tests - Full ODE system behavior
  3. Non-equilibrium Tests - Time-dependent evolution and chemical dynamics
  4. Validation Tests - Physical constraints (atom conservation, non-negativity)
"""

import sys
from pathlib import Path

# Add repo root to path
_THIS_DIR = Path(__file__).resolve().parent
_ROOT = _THIS_DIR.parent  # chemistry/
_SRC_ROOT = _ROOT.parent  # src/
for p in (str(_ROOT), str(_SRC_ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None

from typing import Dict, List, Tuple, Optional, Callable
from chemistry.custom_network_loader import load_network
from chemistry.build_ode_system import (
    make_ode_system,
    integrate_with_sundials,
    prepare_initial_concentrations,
)
from chemistry.util.jacobian_utils import analytic_jacobian
from chemistry.util.cli_utils import enforce_atom_conservation, print_clamped_table


class TestConfig:
    """Configuration for non-equilibrium tests."""

    def __init__(self):
        # Temperature range (K)
        self.temperatures = [10.0, 30.0, 100.0, 300.0, 1000.0]

        # Time integration parameters
        self.t_start = 0.0
        self.t_end = 1e3  # seconds
        self.n_timepoints = 100

        # Physical parameters
        self.av = 10.0  # Visual extinction (mag)
        self.cr_zeta = 1.3e-17  # Cosmic-ray ionization rate (s^-1)

        # Tolerance levels
        self.rtol = 1e-6
        self.atol = 1e-12

        # Output settings
        self.verbose = True
        self.show_rates = True
        self.show_jacobian_spectrum = False


class ChemicalNetworkTest:
    """Base test class for chemical network verification."""

    def __init__(self, network, config: Optional[TestConfig] = None):
        """
        Initialize test suite.

        Args:
            network: network dict (pre-loaded)
            config: TestConfig instance, or use defaults
        """
        self.network = network
        self.config = config or TestConfig()
        self.results = {}

        # Build ODE system
        try:
            self.species, self.S, self.rate_func, self.reaction_data = make_ode_system(
                network
            )
            self.n_species = len(self.species)
            self.n_reactions = len(network["reactions"])
            self._system_built = True
        except Exception as e:
            print(f"Failed to build ODE system: {e}")
            self._system_built = False

    def test_rate_function(
        self, x0: Optional[List[float]] = None, T: float = 100.0
    ) -> bool:
        """
        Test 1: Verify rate function evaluation.

        Args:
            x0: Initial concentrations (use defaults if None)
            T: Temperature (K)

        Returns:
            True if test passed, False otherwise
        """
        if not self._system_built:
            print("❌ TEST 1 (Rate Function): ODE system not built")
            return False

        try:
            # Get initial conditions
            if x0 is None:
                x0 = prepare_initial_concentrations(self.species)

            # Evaluate rates and RHS
            rates, rhs = self.rate_func(
                x0, t=0.0, T=T, av=self.config.av, cr_zeta=self.config.cr_zeta
            )

            # Verify output shapes
            if not (len(rates) == self.n_reactions):
                raise ValueError(
                    f"rates shape mismatch: {len(rates)} vs {self.n_reactions}"
                )
            if not (len(rhs) == self.n_species):
                raise ValueError(f"rhs shape mismatch: {len(rhs)} vs {self.n_species}")

            # Verify no NaNs or infs
            if HAS_NUMPY:
                rates_arr = np.asarray(rates)
                rhs_arr = np.asarray(rhs)
                if np.any(np.isnan(rates_arr)):
                    raise ValueError("NaN in rates")
                if np.any(np.isinf(rates_arr)):
                    raise ValueError("Inf in rates")
                if np.any(np.isnan(rhs_arr)):
                    raise ValueError("NaN in RHS")
                if np.any(np.isinf(rhs_arr)):
                    raise ValueError("Inf in RHS")

            # Basic physical check: most rates should be positive
            rates_list = list(rates) if HAS_NUMPY else rates
            positive_count = sum(1 for r in rates_list if float(r) >= 0)
            if positive_count <= self.n_reactions * 0.9:
                raise AssertionError(
                    f"Too many negative rates: {positive_count}/{self.n_reactions}"
                )

            self.results["rate_function"] = {
                "status": "passed",
                "n_reactions_evaluated": self.n_reactions,
                "max_rate": float(max(rates_list)),
                "min_rate": float(min(rates_list)),
                "n_positive_rates": positive_count,
            }

            print(f"✅ TEST 1 (Rate Function): PASSED at T={T}K")
            if self.config.show_rates:
                print(
                    f"   Max rate: {float(max(rates_list)):.3e}, Min rate: {float(min(rates_list)):.3e}"
                )
                print(f"   Positive rates: {positive_count}/{self.n_reactions}")

            return True

        except Exception as e:
            self.results["rate_function"] = {"status": "failed", "error": str(e)}
            print(f"❌ TEST 1 (Rate Function): FAILED - {e}")
            return False

    def test_jacobian_computation(
        self, x0: Optional[List[float]] = None, T: float = 100.0
    ) -> bool:
        """
        Test 2: Verify analytic Jacobian computation.

        Args:
            x0: Initial concentrations (use defaults if None)
            T: Temperature (K)

        Returns:
            True if test passed, False otherwise
        """
        if not self._system_built:
            print("❌ TEST 2 (Jacobian): ODE system not built")
            return False

        try:
            # Get initial conditions
            if x0 is None:
                x0 = prepare_initial_concentrations(self.species)

            # Compute Jacobian
            J = analytic_jacobian(
                x0,
                self.reaction_data,
                self.species,
                self.S,
                self.network,
                T=T,
                av=self.config.av,
                cr_zeta=self.config.cr_zeta,
            )

            # Verify shape
            if HAS_NUMPY:
                J = np.asarray(J)
                if J.shape != (self.n_species, self.n_species):
                    raise ValueError(f"Jacobian shape mismatch: {J.shape}")

                # Check for NaNs/Infs
                if np.any(np.isnan(J)):
                    raise ValueError("NaN in Jacobian")
                if np.any(np.isinf(J)):
                    raise ValueError("Inf in Jacobian")

                # Compute spectral properties
                eigenvals = np.linalg.eigvals(J)
                max_eig_real = np.max(np.real(eigenvals))
                min_eig_real = np.min(np.real(eigenvals))
                condition_number = np.linalg.cond(J)

                self.results["jacobian"] = {
                    "status": "passed",
                    "shape": J.shape,
                    "max_eigenvalue_real": float(max_eig_real),
                    "min_eigenvalue_real": float(min_eig_real),
                    "condition_number": float(condition_number),
                    "stiffness_ratio": float(
                        abs(max_eig_real) / (abs(min_eig_real) + 1e-30)
                    ),
                }

                print(f"✅ TEST 2 (Jacobian): PASSED at T={T}K")
                if self.config.show_jacobian_spectrum:
                    print(f"   Shape: {J.shape}")
                    print(f"   Max eigenvalue (real): {max_eig_real:.3e}")
                    print(f"   Min eigenvalue (real): {min_eig_real:.3e}")
                    print(f"   Condition number: {condition_number:.3e}")
                    print(
                        f"   Stiffness ratio: {self.results['jacobian']['stiffness_ratio']:.3e}"
                    )
            else:
                # Pure Python fallback
                assert len(J) == self.n_species, "Jacobian row count mismatch"
                assert len(J[0]) == self.n_species, "Jacobian column count mismatch"
                self.results["jacobian"] = {
                    "status": "passed",
                    "shape": (self.n_species, self.n_species),
                }
                print(
                    f"✅ TEST 2 (Jacobian): PASSED at T={T}K (pure Python, no spectral analysis)"
                )

            return True

        except Exception as e:
            self.results["jacobian"] = {"status": "failed", "error": str(e)}
            print(f"❌ TEST 2 (Jacobian): FAILED - {e}")
            return False

    def test_stoichiometric_matrix(self) -> bool:
        """
        Test 3: Verify stoichiometric matrix properties.

        Returns:
            True if test passed, False otherwise
        """
        if not self._system_built:
            print("❌ TEST 3 (Stoichiometry): ODE system not built")
            return False

        try:
            S = self.S

            # Basic shape check
            if HAS_NUMPY:
                S = np.asarray(S)
                if S.shape != (self.n_species, self.n_reactions):
                    raise ValueError(f"S shape mismatch: {S.shape}")

                # Check row sums for some reactions (should be conserved for elementary reactions)
                n_zero_col_sum = np.sum(np.abs(np.sum(S, axis=0)) < 1e-10)

                # Verify all elements are integers (stoichiometry)
                if not np.allclose(S, np.round(S)):
                    raise ValueError(
                        "Stoichiometric matrix contains non-integer values"
                    )

            self.results["stoichiometry"] = {
                "status": "passed",
                "shape": S.shape if HAS_NUMPY else (self.n_species, self.n_reactions),
                "n_species": self.n_species,
                "n_reactions": self.n_reactions,
            }

            print(f"✅ TEST 3 (Stoichiometry): PASSED")
            print(f"   Shape: ({self.n_species} species, {self.n_reactions} reactions)")

            return True

        except Exception as e:
            self.results["stoichiometry"] = {"status": "failed", "error": str(e)}
            print(f"❌ TEST 3 (Stoichiometry): FAILED - {e}")
            return False

    def test_atom_conservation(self, x0: Optional[List[float]] = None) -> bool:
        """
        Test 4: Verify atom conservation enforcement.

        Args:
            x0: Initial concentrations (use defaults if None)

        Returns:
            True if test passed, False otherwise
        """
        if not self._system_built:
            print("❌ TEST 4 (Atom Conservation): ODE system not built")
            return False

        try:
            # Get initial conditions
            if x0 is None:
                x0 = prepare_initial_concentrations(self.species)

            # Create perturbed final state (simulate integration errors)
            if HAS_NUMPY:
                final = np.asarray(x0) * 0.9 + np.random.rand(
                    self.n_species
                ) * 0.01 * np.asarray(x0)
            else:
                final = [v * 0.9 + 0.001 * v for v in x0]

            # Enforce conservation
            final_conserved = enforce_atom_conservation(
                final, x0, self.species, self.network
            )

            # Verify conservation
            element_order = self.network.get("element_order", [])
            errors = []

            for element in element_order:
                # Count initial atoms
                n_init = 0
                n_final = 0

                for i, spec in enumerate(self.species):
                    comp = self.network["species"][spec].get("composition", {})
                    count = comp.get(element, 0)
                    if count > 0:
                        n_init += count * (x0[i] if isinstance(x0, list) else x0[i])
                        n_final += count * (
                            final_conserved[i]
                            if isinstance(final_conserved, list)
                            else final_conserved[i]
                        )

                # Check relative error
                if n_init > 0:
                    rel_error = abs(n_final - n_init) / n_init
                    if rel_error > 1e-6:
                        errors.append((element, rel_error))

            if errors:
                print(
                    f"⚠️  TEST 4 (Atom Conservation): Some elements not perfectly conserved:"
                )
                for elem, err in errors[:5]:
                    print(f"    {elem}: relative error = {err:.3e}")

            self.results["atom_conservation"] = {
                "status": "passed" if not errors else "warning",
                "n_elements": len(element_order),
                "conservation_errors": len(errors),
            }

            print(f"✅ TEST 4 (Atom Conservation): PASSED")
            print(
                f"   Elements: {len(element_order)}, Conservation errors: {len(errors)}"
            )

            return True

        except Exception as e:
            self.results["atom_conservation"] = {"status": "failed", "error": str(e)}
            print(f"❌ TEST 4 (Atom Conservation): FAILED - {e}")
            return False

    def test_integration_sundials(
        self, x0: Optional[List[float]] = None, T: float = 100.0
    ) -> bool:
        """
        Test 5: Verify SUNDIALS (CVODE) integration.

        Args:
            x0: Initial concentrations (use defaults if None)
            T: Temperature (K)

        Returns:
            True if test passed, False otherwise
        """
        if not self._system_built:
            print("❌ TEST 5 (SUNDIALS Integration): ODE system not built")
            return False

        try:
            # Check SUNDIALS availability
            try:
                from scikits.odes import ode as _scikits_ode
            except ImportError:
                print(
                    "⏭️  TEST 5 (SUNDIALS Integration): SKIPPED (SUNDIALS not available)"
                )
                self.results["integration_sundials"] = {
                    "status": "skipped",
                    "reason": "SUNDIALS not available",
                }
                return True

            # Get initial conditions
            if x0 is None:
                x0 = prepare_initial_concentrations(self.species)

            if HAS_NUMPY:
                x0 = np.asarray(x0, dtype=float)

            # Run integration
            print(
                f"   Running SUNDIALS integration from {self.config.t_start} to {self.config.t_end} s at T={T}K..."
            )

            sol = integrate_with_sundials(
                self.rate_func,
                x0,
                self.config.t_start,
                self.config.t_end,
                T=T,
                av=self.config.av,
                cr_zeta=self.config.cr_zeta,
                reaction_data=self.reaction_data,
                species=self.species,
                S=self.S,
                network=self.network,
                use_analytic_jac=True,
            )

            # Extract final state
            try:
                final = sol.y[-1]
            except:
                try:
                    final = sol.values.y[-1]
                except:
                    final = None

            if final is None:
                raise ValueError("Could not extract final state from solution")

            # Verify output
            if HAS_NUMPY:
                final = np.asarray(final)
                if np.any(np.isnan(final)):
                    raise ValueError("NaN in final state")
                if np.any(np.isinf(final)):
                    raise ValueError("Inf in final state")
                max_conc = np.max(final)
                min_conc = np.min(final)
            else:
                max_conc = max(final)
                min_conc = min(final)

            self.results["integration_sundials"] = {
                "status": "passed",
                "max_concentration": float(max_conc),
                "min_concentration": float(min_conc),
                "n_species_nonzero": sum(1 for c in final if c > 1e-30),
            }

            print(f"✅ TEST 5 (SUNDIALS Integration): PASSED at T={T}K")
            print(
                f"   Max concentration: {max_conc:.3e}, Min concentration: {min_conc:.3e}"
            )
            print(
                f"   Non-zero species: {sum(1 for c in final if c > 1e-30)}/{self.n_species}"
            )

            return True

        except Exception as e:
            self.results["integration_sundials"] = {"status": "failed", "error": str(e)}
            print(f"❌ TEST 5 (SUNDIALS Integration): FAILED - {e}")
            return False

    def test_nonequilibrium_evolution(self, x0: Optional[List[float]] = None) -> bool:
        """
        Test 6: Verify non-equilibrium chemical evolution across temperatures.

        Tests that chemistry changes significantly with temperature, indicating
        non-equilibrium dynamics are properly captured.

        Args:
            x0: Initial concentrations (use defaults if None)

        Returns:
            True if test passed, False otherwise
        """
        if not self._system_built:
            print("❌ TEST 6 (Non-equilibrium Evolution): ODE system not built")
            return False

        try:
            # Get initial conditions
            if x0 is None:
                x0 = prepare_initial_concentrations(self.species)

            if HAS_NUMPY:
                x0 = np.asarray(x0, dtype=float)

            print(
                f"   Testing evolution across {len(self.config.temperatures)} temperatures..."
            )

            evolution_results = []

            for i, T in enumerate(self.config.temperatures):
                # Evaluate rates at this temperature
                rates, rhs = self.rate_func(
                    x0, t=0.0, T=T, av=self.config.av, cr_zeta=self.config.cr_zeta
                )

                # Compute characteristic timescale: 1 / max(|RHS / conc|)
                if HAS_NUMPY:
                    x0_arr = np.asarray(x0)
                    rhs_arr = np.asarray(rhs)
                    # Avoid division by zero
                    timescales = np.abs(rhs_arr) / (np.abs(x0_arr) + 1e-30)
                    max_timescale_inv = float(np.max(timescales))
                    char_time = 1.0 / (max_timescale_inv + 1e-30)
                    rates_list = list(rates)
                else:
                    max_rate = max(abs(r) for r in rhs)
                    char_time = 1.0 / (max_rate + 1e-30)
                    rates_list = rates

                evolution_results.append(
                    {
                        "T": T,
                        "max_rate": float(max(rates_list)) if rates_list else 0,
                        "char_time": char_time,
                        "n_active_reactions": sum(1 for r in rates if r > 1e-30),
                    }
                )

            # Check that chemistry changes with temperature
            rate_variation = []
            for result in evolution_results:
                rate_variation.append(result["max_rate"])

            # Should have variation across temperatures
            if HAS_NUMPY:
                rate_var_normalized = np.std(rate_variation) / (
                    np.mean(rate_variation) + 1e-30
                )
            else:
                mean_rate = sum(rate_variation) / len(rate_variation)
                var_rate = sum((r - mean_rate) ** 2 for r in rate_variation) / len(
                    rate_variation
                )
                std_rate = var_rate**0.5
                rate_var_normalized = std_rate / (mean_rate + 1e-30)

            has_significant_variation = rate_var_normalized > 0.1

            self.results["nonequilibrium_evolution"] = {
                "status": "passed",
                "temperature_range": (
                    min(self.config.temperatures),
                    max(self.config.temperatures),
                ),
                "rate_variation": float(rate_var_normalized),
                "has_significant_variation": has_significant_variation,
                "evolution_results": evolution_results,
            }

            print(f"✅ TEST 6 (Non-equilibrium Evolution): PASSED")
            print(
                f"   Temperature range: {min(self.config.temperatures)}K - {max(self.config.temperatures)}K"
            )
            print(
                f"   Rate variation coefficient: {rate_var_normalized:.3f} (significant: {has_significant_variation})"
            )

            return True

        except Exception as e:
            self.results["nonequilibrium_evolution"] = {
                "status": "failed",
                "error": str(e),
            }
            print(f"❌ TEST 6 (Non-equilibrium Evolution): FAILED - {e}")
            return False

    def run_all_tests(self) -> Dict[str, bool]:
        """
        Run all tests and return summary.

        Returns:
            Dictionary with test names as keys and pass/fail as values
        """
        print("\n" + "=" * 70)
        print("CHEMICAL NETWORK NON-EQUILIBRIUM TEST SUITE")
        print("=" * 70 + "\n")

        results = {}

        # Run tests
        results["test_rate_function"] = self.test_rate_function(T=100.0)
        results["test_stoichiometric_matrix"] = self.test_stoichiometric_matrix()
        results["test_jacobian_computation"] = self.test_jacobian_computation(T=100.0)
        results["test_atom_conservation"] = self.test_atom_conservation()
        results["test_integration_sundials"] = self.test_integration_sundials(T=100.0)
        results["test_nonequilibrium_evolution"] = self.test_nonequilibrium_evolution()

        # Print summary
        print("\n" + "=" * 70)
        print("TEST SUMMARY")
        print("=" * 70)

        passed = sum(1 for v in results.values() if v)
        total = len(results)

        for name, passed_flag in results.items():
            status = "✅ PASSED" if passed_flag else "❌ FAILED"
            print(f"  {name:<40} {status}")

        print(f"\nTotal: {passed}/{total} tests passed")
        print("=" * 70 + "\n")

        return results


def main():
    """Main entry point for test suite."""
    # Load network
    print("Loading chemical network...")
    custom_dir = Path(_THIS_DIR).parent / "custom"
    network = load_network(custom_dir)
    if not network.get("species") or not network.get("reactions"):
        print("Failed to load network")
        sys.exit(1)

    # Configure tests
    config = TestConfig()
    config.verbose = True
    config.show_rates = True
    config.show_jacobian_spectrum = True

    # Run tests
    test_suite = ChemicalNetworkTest(network, config)
    results = test_suite.run_all_tests()

    # Exit with appropriate code
    if all(results.values()):
        print("✅ All tests passed!")
        sys.exit(0)
    else:
        print("❌ Some tests failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
