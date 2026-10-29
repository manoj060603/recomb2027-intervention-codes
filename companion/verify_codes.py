#!/usr/bin/env python3
"""Exact finite intervention-code certificates; no data or simulation.

Only Python's standard library is used. Vertex labels in JSON are 1-based.
The rank field is Q, which gives exactly the same rank as R for these
integer matrices. Decimal analytic constants are explanatory evaluations;
the symbolic formulas in Appendix K are the authoritative bounds.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
import json
import math
from pathlib import Path
import sys
import unittest


def rref(matrix):
    """Return exact RREF and pivot columns without modifying the input."""
    a = [[Fraction(x) for x in row] for row in matrix]
    if not a:
        return a, []
    width = len(a[0])
    if any(len(row) != width for row in a):
        raise ValueError("Ragged matrix")
    pivot_columns = []
    next_row = 0
    for col in range(width):
        pivot = next((i for i in range(next_row, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[next_row], a[pivot] = a[pivot], a[next_row]
        divisor = a[next_row][col]
        a[next_row] = [x / divisor for x in a[next_row]]
        for row in range(len(a)):
            if row != next_row and a[row][col]:
                multiplier = a[row][col]
                a[row] = [x - multiplier * y for x, y in zip(a[row], a[next_row])]
        pivot_columns.append(col)
        next_row += 1
        if next_row == len(a):
            break
    return a, pivot_columns


def determinant(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    n = len(a)
    if any(len(row) != n for row in a):
        raise ValueError("Determinant requires a square matrix")
    result = Fraction(1)
    for col in range(n):
        pivot = next((i for i in range(col, n) if a[i][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            result = -result
        diagonal = a[col][col]
        result *= diagonal
        for row in range(col + 1, n):
            factor = a[row][col] / diagonal
            for j in range(col + 1, n):
                a[row][j] -= factor * a[col][j]
            a[row][col] = Fraction(0)
    return result


def primitive_null_vector(matrix):
    reduced, pivots = rref(matrix)
    width = len(matrix[0])
    free = next((j for j in range(width) if j not in pivots), None)
    if free is None:
        raise ValueError("Matrix has no nonzero null vector")
    v = [Fraction(0) for _ in range(width)]
    v[free] = Fraction(1)
    for row, pivot in enumerate(pivots):
        v[pivot] = -reduced[row][free]
    common = math.lcm(*(x.denominator for x in v))
    integers = [int(x * common) for x in v]
    divisor = math.gcd(*integers)
    integers = [x // divisor for x in integers]
    first = next(x for x in integers if x)
    if first < 0:
        integers = [-x for x in integers]
    if any(sum(Fraction(a) * b for a, b in zip(row, integers)) for row in matrix):
        raise AssertionError("Internal null-vector certificate failure")
    return integers


def exact_string(x):
    x = Fraction(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def validate_problem(problem):
    p = problem.get("p")
    r = problem.get("r")
    if isinstance(p, bool) or not isinstance(p, int) or p < 1:
        raise ValueError("p must be a positive integer")
    if isinstance(r, bool) or not isinstance(r, int) or not 0 <= r < p:
        raise ValueError("r must be an integer in {0,...,p-1}")
    edges = set()
    for pair in problem.get("reciprocal_edges", []):
        if len(pair) != 2 or any(isinstance(x, bool) or not isinstance(x, int) for x in pair):
            raise ValueError("Each reciprocal edge must contain two integer labels")
        u, v = pair
        if not (1 <= u <= p and 1 <= v <= p) or u == v:
            raise ValueError("Reciprocal edges must join distinct vertices in [p]")
        edge = tuple(sorted((u, v)))
        if edge in edges:
            raise ValueError("Duplicate reciprocal edge")
        edges.add(edge)
    target_sets = []
    for targets in problem.get("target_sets", []):
        if any(isinstance(i, bool) or not isinstance(i, int) or not 1 <= i <= p for i in targets):
            raise ValueError("Target labels must be integers in [p]")
        if len(set(targets)) != len(targets):
            raise ValueError("Duplicate label in a target set")
        target_sets.append(set(targets))
    budget = problem.get("target_budget", p)
    if isinstance(budget, bool) or not isinstance(budget, int) or not 1 <= budget <= p:
        raise ValueError("target_budget must be an integer in [1,p]")
    if any(len(targets) > budget for targets in target_sets):
        raise ValueError("A target set exceeds target_budget")
    return p, r, edges, target_sets, budget


def verify_problem(problem):
    p, r, edges, target_sets, budget = validate_problem(problem)
    matrix = [[1] * p] + [[int(i not in targets) for i in range(1, p + 1)]
                               for targets in target_sets]
    report = {
        "name": problem.get("name", "unnamed"),
        "arithmetic": "exact rational; identical rank over the real numbers",
        "p": p,
        "r": r,
        "m": len(target_sets),
        "total_contexts_including_baseline": len(matrix),
        "target_budget": budget,
        "maximum_target_size": max((len(s) for s in target_sets), default=0),
        "targeting_incidences": sum(len(s) for s in target_sets),
        "retention_matrix": matrix,
        "target_codes": [[int(i in targets) for targets in target_sets]
                         for i in range(1, p + 1)],
        "certificates": [],
    }
    for size in range(1, r + 2):
        for clique in combinations(range(1, p + 1), size):
            if not all((u, v) in edges for u, v in combinations(clique, 2)):
                continue
            submatrix = [[row[i - 1] for i in clique] for row in matrix]
            _, pivot_columns = rref(submatrix)
            rank = len(pivot_columns)
            if rank < size:
                report.update({
                    "success": False,
                    "feasible_cliques_checked": len(report["certificates"]) + 1,
                    "deficient_clique": list(clique),
                    "deficient_clique_size": size,
                    "deficient_rank": rank,
                    "null_vector": primitive_null_vector(submatrix),
                    "null_vector_coordinate_order": list(clique),
                    "minimal_obstruction": True,
                })
                return report
            transpose = [list(column) for column in zip(*submatrix)]
            _, selected_rows = rref(transpose)
            selected_rows = selected_rows[:size]
            minor = [[submatrix[i][j] for j in range(size)] for i in selected_rows]
            det = determinant(minor)
            if not det:
                raise AssertionError("Internal success-minor certificate failure")
            report["certificates"].append({
                "clique": list(clique),
                "rank": rank,
                "minor_context_rows": selected_rows,
                "minor_determinant": exact_string(det),
            })
    report.update({"success": True,
                   "feasible_cliques_checked": len(report["certificates"])})
    return report


def complete_edges(p):
    return [list(pair) for pair in combinations(range(1, p + 1), 2)]


def analytic_constants():
    c = 9 * math.sqrt(3)
    d = 27 * math.sqrt(3) / 2
    q = 729 / 2
    a_roots = [(7 + sign * math.sqrt(33)) / 4 for sign in (-1, 1)]
    b_roots = [(5 + sign * math.sqrt(17)) / 4 for sign in (-1, 1)]
    a_star = max(2 * math.sqrt(u) * abs(u - 2) * math.exp(-u) for u in a_roots)
    b_star = max(2 * math.sqrt(u) * abs(1 - u) * math.exp(-u) for u in b_roots)

    def radial_maximum(center, exponent):
        linear = center + 2 * exponent + 2
        discriminant = linear ** 2 - 8 * exponent * center
        candidates = [(linear + sign * math.sqrt(discriminant)) / 2
                      for sign in (-1, 1)]
        return max(abs(center - s) * s ** exponent * math.exp(-s / 2)
                   for s in candidates)

    bpsi = c * math.exp(-1.5)
    dpsi = c * math.sqrt(math.exp(-4) + 3 * math.exp(-2))
    bq = q * math.exp(-3)
    lq = q * math.sqrt(a_star ** 2 * math.exp(-6) + 3 * b_star ** 2 * math.exp(-4))
    br = d / 8 * radial_maximum(8, 2)
    lr = 2 * d * (radial_maximum(8, 1.5) / 4 + radial_maximum(10, 2.5) / 16)
    kbar = math.sqrt(10) * lq
    sigma = math.sqrt(3) - 1
    eps = sigma / (2 * kbar)
    planning = []
    for moment_count in (14, 20):
        logarithm = math.log(2 * moment_count / 0.05)
        cm = math.sqrt(10) * bq * math.sqrt(2 * logarithm)
        cb = 2 * br * math.sqrt(2 * logarithm)
        threshold = (2 * cm / sigma) ** 2
        planning.append({
            "tagged_moment_count": moment_count,
            "confidence_failure_probability": 0.05,
            "matrix_allowance_prefactor": cm,
            "vector_allowance_prefactor": cb,
            "reference_strict_sample_threshold_per_context": threshold,
            "reference_ideal_minimum_integer_per_context": math.floor(threshold) + 1,
            "reference_ideal_total_cells_four_contexts": 4 * (math.floor(threshold) + 1),
            "sufficient_ball_strict_sample_threshold_per_context": 4 * threshold,
            "sufficient_ball_ideal_minimum_integer_per_context": math.floor(4 * threshold) + 1,
            "sufficient_ball_ideal_total_cells_four_contexts": 4 * (math.floor(4 * threshold) + 1),
        })
    return {
        "purpose": "deterministic evaluation of the analytic formulas in Appendix K",
        "decimal_values_are_approximations_not_certified_intervals": True,
        "authoritative_bounds": "symbolic formulas in Appendix K",
        "exact_identities": {
            "gram_matrix": "I_4",
            "F_squared_Gaussian_L2_norm": "4/3",
            "good_reference_matrix": [[1, 1, 1, 1], [1, 0, 0, 1],
                                      [1, 0, 1, 0], [1, 1, 0, 0]],
            "good_reference_singular_values": ["sqrt(3)-1", "1", "1", "sqrt(3)+1"],
            "bad_null_vector": [1, -1, -1, 1],
            "equal_cost_bad_reference_matrix": [[1, 1, 1, 1], [1, 0, 1, 0],
                                                [1, 1, 0, 0], [0, 1, 0, 1]],
            "equal_cost_bad_null_vector": [1, -1, -1, 1],
            "equal_cost_bad_rank": 3,
            "nonzero_tagged_moments_for_good_and_equal_cost_bad": 14,
        },
        "symbolic_definitions": {
            "D": "27*sqrt(3)/2",
            "Q": "729/2",
            "B_psi": "9*sqrt(3)*exp(-3/2)",
            "D_psi": "9*sqrt(3)*sqrt(exp(-4)+3*exp(-2))",
            "A_star": "max_{u=(7+-sqrt(33))/4} 2*sqrt(u)*abs(u-2)*exp(-u)",
            "B_star": "max_{u=(5+-sqrt(17))/4} 2*sqrt(u)*abs(1-u)*exp(-u)",
            "R(c,a)": "max_{s=(c+2*a+2+-sqrt((c+2*a+2)^2-8*a*c))/2} abs(c-s)*s^a*exp(-s/2)",
            "B_q": "(729/2)*exp(-3)",
            "L_q": "Q*sqrt(A_star^2*exp(-6)+3*B_star^2*exp(-4))",
            "B_r": "(D/8)*R(8,2)",
            "L_r": "2*D*(R(8,3/2)/4+R(10,5/2)/16)",
            "K_bar": "sqrt(10)*L_q",
            "epsilon_ex": "(sqrt(3)-1)/(2*K_bar)",
            "rho_ex": "epsilon_ex/B_psi",
            "optional_contraction_coefficient_radius": "1/(2*D_psi)",
        },
        "bound_interpretations": {
            "B_r": "exact supremum norm of r for the declared observable",
            "L_q_and_L_r": "proved upper bounds, not exact multivariate Lipschitz constants",
            "K_bar": "conservative upper bound for the family matrix perturbation constant",
            "epsilon_ex_and_rho_ex": "sufficient local radii, not identifiability thresholds",
        },
        "approximate_values": {
            "sigma_star": sigma, "B_psi": bpsi, "D_psi": dpsi,
            "A_star": a_star, "B_star": b_star,
            "B_q": bq, "L_q": lq, "B_r": br, "L_r": lr,
            "K_bar": kbar, "epsilon_ex": eps, "rho_ex": eps / bpsi,
            "optional_contraction_coefficient_radius": 1 / (2 * dpsi),
        },
        "costs": [
            {"design": "bad", "m": 2, "maximum_target_size": 2, "targeting_incidences": 4},
            {"design": "equal_cost_bad", "m": 3, "maximum_target_size": 2, "targeting_incidences": 6},
            {"design": "good", "m": 3, "maximum_target_size": 2, "targeting_incidences": 6},
            {"design": "isolation", "m": 3, "maximum_target_size": 3, "targeting_incidences": 9},
        ],
        "idealized_equal_count_planning": {
            "assumptions": [
                "fixed tests and equal n iid stationary joint cells in each of four used contexts",
                "no additional drift, sampling-law, evaluation or matrix-entry error allowances",
                "ideal exact lower singular-value certificate for the displayed minimum counts",
                "14 counts only nonzero tagged moments; 20 is the generic N^2+N upper count",
            ],
            "symbolic_formulas": {
                "ell": "log(2*L/alpha)",
                "c_M": "sqrt(10)*B_q*sqrt(2*ell)",
                "c_b": "2*B_r*sqrt(2*ell)",
                "beta_M": "c_M/sqrt(n)",
                "beta_b": "c_b/sqrt(n)",
                "reference_strict_threshold": "n > (2*c_M/sigma_star)^2",
                "sufficient_ball_strict_threshold": "n > 4*(2*c_M/sigma_star)^2",
            },
            "numerical_tightness_requirement": "s(Mhat)-eta_s <= shat <= s(Mhat); sufficient reference passage condition: 2*beta_M+eta_s < sigma_star",
            "qualifications": [
                "floating-point evaluations are not certified planning intervals",
                "minimum integers have very little slack; account for numerical tightness explicitly",
                "sufficient planning counts are not necessary, minimax, observed passage predictions, or useful-accuracy guarantees",
            ],
            "approximate_values": planning,
        },
    }


class ExactCertificateTests(unittest.TestCase):
    def problem(self, p, r, target_sets, edges=None, budget=None):
        return {"p": p, "r": r, "target_sets": target_sets,
                "reciprocal_edges": complete_edges(p) if edges is None else edges,
                "target_budget": p if budget is None else budget}

    def test_bad_four_gene_minimal_obstruction(self):
        result = verify_problem(self.problem(4, 3, [[2, 4], [3, 4]], budget=2))
        self.assertFalse(result["success"])
        self.assertEqual(result["deficient_clique"], [1, 2, 3, 4])
        self.assertEqual(result["null_vector"], [1, -1, -1, 1])
        self.assertEqual(result["deficient_rank"], 3)

    def test_good_four_gene_all_cliques(self):
        result = verify_problem(self.problem(4, 3, [[2, 3], [2, 4], [3, 4]], budget=2))
        self.assertTrue(result["success"])
        self.assertEqual(result["feasible_cliques_checked"], 15)
        self.assertEqual(abs(int(result["certificates"][-1]["minor_determinant"])), 2)

    def test_equal_cost_bad_four_gene_exact_collision(self):
        bad = verify_problem(self.problem(4, 3, [[2, 4], [3, 4], [1, 3]], budget=2))
        good = verify_problem(self.problem(4, 3, [[2, 3], [2, 4], [3, 4]], budget=2))
        self.assertFalse(bad["success"])
        self.assertEqual(bad["deficient_clique"], [1, 2, 3, 4])
        self.assertEqual(bad["deficient_rank"], 3)
        self.assertEqual(bad["null_vector"], [1, -1, -1, 1])
        self.assertEqual(bad["retention_matrix"],
                         [[1, 1, 1, 1], [1, 0, 1, 0], [1, 1, 0, 0], [0, 1, 0, 1]])
        self.assertEqual(determinant(bad["retention_matrix"]), 0)
        for field in ("m", "total_contexts_including_baseline", "maximum_target_size", "targeting_incidences"):
            self.assertEqual(bad[field], good[field])
        for report in (bad, good):
            # Four tagged r moments plus one q moment per retained component.
            self.assertEqual(4 + sum(map(sum, report["retention_matrix"])), 14)

    def test_six_gene_order_three(self):
        result = verify_problem(self.problem(6, 3, [[2, 6], [3, 6], [4, 6], [5, 6]], budget=2))
        self.assertTrue(result["success"])
        self.assertEqual(result["feasible_cliques_checked"], 56)

    def test_six_gene_full_order_detects_real_circuit(self):
        result = verify_problem(self.problem(6, 5, [[2, 6], [3, 6], [4, 6], [5, 6]], budget=2))
        self.assertFalse(result["success"])
        self.assertEqual(result["deficient_clique"], [1, 2, 3, 4, 5, 6])
        self.assertEqual(result["null_vector"], [3, -1, -1, -1, -1, 1])

    def test_one_gene_zero_order(self):
        result = verify_problem(self.problem(1, 0, [], edges=[], budget=1))
        self.assertTrue(result["success"])
        self.assertEqual(result["feasible_cliques_checked"], 1)

    def test_zero_order_requires_only_singletons(self):
        result = verify_problem(self.problem(4, 0, [], budget=4))
        self.assertTrue(result["success"])
        self.assertEqual(result["feasible_cliques_checked"], 4)

    def test_edgeless_full_order_requires_only_singletons(self):
        result = verify_problem(self.problem(4, 3, [], edges=[], budget=4))
        self.assertTrue(result["success"])
        self.assertEqual(result["feasible_cliques_checked"], 4)

    def test_singleton_target_budget(self):
        self.assertTrue(verify_problem(self.problem(4, 3, [[1], [2], [3]], budget=1))["success"])

    def test_full_target_budget_does_not_imply_rank(self):
        result = verify_problem(self.problem(4, 3, [[1, 2, 3, 4]], budget=4))
        self.assertFalse(result["success"])
        self.assertEqual(result["deficient_clique"], [1, 2])
        self.assertEqual(result["null_vector"], [1, -1])

    def test_wide_matrix_rank_and_null_vector(self):
        matrix = [[1, 1, 1, 1]]
        self.assertEqual(len(rref(matrix)[1]), 1)
        self.assertEqual(primitive_null_vector(matrix), [1, -1, 0, 0])

    def test_budget_violation_rejected(self):
        with self.assertRaises(ValueError):
            verify_problem(self.problem(4, 3, [[2, 3]], budget=1))

    def test_good_gram_characteristic_polynomial(self):
        a = [[1, 1, 1, 1], [1, 0, 0, 1], [1, 0, 1, 0], [1, 1, 0, 0]]
        gram = [[sum(row[i] * row[j] for row in a) for j in range(4)] for i in range(4)]
        self.assertEqual(gram, [[4, 2, 2, 2], [2, 2, 1, 1], [2, 1, 2, 1], [2, 1, 1, 2]])
        # Five exact evaluations verify (lambda-1)^2*(lambda^2-8lambda+4),
        # since the two polynomials are monic degree four and five points suffice.
        for lam in range(5):
            shifted = [[gram[i][j] - (lam if i == j else 0) for j in range(4)] for i in range(4)]
            self.assertEqual(determinant(shifted), (lam - 1) ** 2 * (lam ** 2 - 8 * lam + 4))


def write_json(value, destination=None):
    rendered = json.dumps(value, indent=2, sort_keys=True) + "\n"
    if destination is None:
        sys.stdout.write(rendered)
    else:
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(rendered, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--input", type=Path, help="JSON problem to verify")
    mode.add_argument("--all-examples", action="store_true")
    mode.add_argument("--constants", action="store_true")
    mode.add_argument("--self-test", action="store_true")
    parser.add_argument("--output", type=Path, help="JSON output file for --input or --constants")
    parser.add_argument("--output-dir", type=Path, help="certificate directory for --all-examples")
    args = parser.parse_args()
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(ExactCertificateTests)
        result = unittest.TextTestRunner(stream=sys.stderr, verbosity=2).run(suite)
        write_json({"deterministic_self_checks": result.testsRun,
                    "success": result.wasSuccessful(),
                    "uses_randomness_or_observations": False})
        return 0 if result.wasSuccessful() else 1
    if args.constants:
        write_json(analytic_constants(), args.output)
        return 0
    if args.input:
        try:
            problem = json.loads(args.input.read_text(encoding="utf-8"))
            report = verify_problem(problem)
        except (OSError, ValueError, TypeError, KeyError) as exc:
            parser.error(str(exc))
        write_json(report, args.output)
        return 0
    base = Path(__file__).resolve().parent
    reports = {}
    for source in sorted((base / "examples").glob("*.json")):
        report = verify_problem(json.loads(source.read_text(encoding="utf-8")))
        reports[source.stem] = report
        if args.output_dir:
            write_json(report, args.output_dir / source.name)
    if args.output_dir is None:
        write_json(reports)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
