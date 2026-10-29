# Deterministic algebra companion

This companion verifies the finite intervention-code criterion and
evaluates the explicit analytic constants in Appendix K. It contains no
datasets, simulated trajectories, randomized computations, experimental
performance claims, or random seeds. No third-party packages are required.

The only runtime requirement is Python 3.9 or later. All rank, determinant,
and null-vector operations use `fractions.Fraction`. The integer retention
matrices therefore have exactly their real rank, rather than a rank inferred
from a numerical tolerance or from arithmetic over a finite field.

From the paper project directory, reproduce every saved certificate with:

```sh
python3 companion/verify_codes.py --self-test
python3 companion/verify_codes.py --all-examples --output-dir companion/certificates
python3 companion/verify_codes.py --constants --output companion/certificates/analytic_constants.json
```

Verify a single declared graph and design with:

```sh
python3 companion/verify_codes.py --input companion/examples/k4_good.json
```

The example JSON files declare the reciprocal graph, interaction order,
target sets, and maximum target-set budget. Gene labels start at one; context
row zero in a certificate is baseline. The verifier enumerates all nonempty
reciprocal cliques of size at most `r + 1`, in increasing size and then
lexicographic order. A success certificate lists a nonzero square minor for
each clique. A failure certificate supplies the first deficient clique and
an integer null vector, aligned with its listed vertices. Its proper
subsets were already checked, so this first obstruction is minimal.

The examples are the invalid four-gene square codes, the valid four-gene
tetrahedron codes, and the six-gene four-bit construction. The self-checks
also cover `p = 1`, zero interaction order, an edgeless reciprocal graph,
wide matrices, singleton targeting, the full targeting budget, budget
violations, and failure of the six-gene construction when the interaction
order is increased to five.

The companion does not enumerate optimal designs or prove universal
complexity bounds. Explicit clique enumeration is finite and can be
exponential in the number of genes; no efficient or uniform runtime claim
is made. The separate countable test-selection procedure in Appendix F is
conditional on effective integration and evaluation contracts and has no
uniform stopping bound. Failure of a finite test pool is inconclusive.

`analytic_constants.json` reproduces the symbolic formulas and supplies
deterministic floating-point evaluations for orientation. Those decimal
values are not certified intervals and must not be substituted for
conservative numerical error bounds. The exact formulas in Appendix K
define the proved envelopes and sufficient radius. The analytic basis is

```text
psi_j(x) = 9 sqrt(3) e_j product_{k != j} x_k exp(-|x|^2/2)
phi(x) = (27 sqrt(3)/2) product_k x_k exp(-|x|^2/2).
```

Reference moments, orthonormality, the common divergence, and the singular
values are proved by Gaussian integration and exact linear algebra in
Appendix K. They do not require numerical quadrature. For other basis/test
families, the certified integration and Gram-matrix conditions of Appendix F
must be supplied; arbitrary black-box smooth functions do not satisfy an
effective contract merely by being smooth.

Finite-sample use is separate from these deterministic outputs. The proved
statistical certificate requires the stated independent stationary joint
cells within each used context, fixed tests or valid selection control,
per-context sample counts, certified moment/evaluation errors, and an
observed positive singular-value margin. Reusing a cell across fixed
moments is permitted. Trajectory dependence, nonstationarity, measurement
distortion, misspecification, and numerical bias require their own stated
assumptions or allowances. No sampled certificate is supplied here.

Public companion repository: https://github.com/manoj060603/recomb2027-intervention-codes

This directory is the reproducibility package for the displayed algebraic
checks. It contains no empirical study. The exact source snapshot is recorded
in the submission package's audit/release-manifest.json.
