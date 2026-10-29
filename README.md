# Intervention codes for nonlinear regulatory diffusions

Deterministic companion to **Intervention Codes and Local Identifiability of
Nonlinear Regulatory Diffusions**, by Manoj Saravanan, Department of ECE,
Virginia Tech (RECOMB 2027 submission manuscript).

This repository verifies finite designs and reproduces analytic formulas.
It contains no data, simulations, fitted biological network, or empirical
validation. The manuscript supplies the mathematical proofs.

## Reproduce the exact certificates

Python 3.9 or later; standard library only. Run from the repository root:

```sh
python3 companion/verify_codes.py --self-test
python3 companion/verify_codes.py --all-examples --output-dir companion/certificates
python3 companion/verify_codes.py --constants --output companion/certificates/analytic_constants.json
```

The first command performs 13 deterministic self-checks. The next commands
reproduce the original invalid four-gene square, the valid four-gene design,
an equal-cost invalid four-gene design, the six-gene order-three design, and
the sharper analytic constants from Appendix K.

All ranks and null vectors use exact rational arithmetic. Success certificates
list nonzero minors for every feasible reciprocal clique. A failure certificate
lists a deficient clique and an exact null vector. Decimal evaluations of
analytic constants are explanatory approximations, not certified intervals.
The improved `B_r` is the exact observable supremum norm; `L_q` and `L_r`
are Lipschitz upper bounds, and `K_bar` is a conservative family bound.
Their sufficient radii are not identification thresholds. The equal-cost
positive/negative comparison is exact algebra, without sampled observations.
The constants JSON distinguishes fourteen nonzero tagged moments from the
generic twenty-moment count in an explicitly idealized planning illustration;
its counts are not necessary sample sizes or useful-accuracy guarantees.

See [companion/README.md](companion/README.md) for the input contract, graph
and target conventions, fixed-family statistical assumptions, integration
qualifications, and scope. No optimal-design runtime or uniform sample bound
is claimed. No seeds are needed because these computations use no randomness.

## Reproduce the schematic code diagram

The optional figure generator requires ReportLab, tested with the version
recorded in the manuscript's source package. From the repository root:

```sh
python3 -m pip install reportlab
python3 figures/generate_code_geometry.py
```

It writes `figures/code-geometry.pdf`, a vector diagram of declared codewords;
it does not plot observations or simulated results.

`MANIFEST.sha256` records the published implementation, inputs, and certificates.
Repository: https://github.com/manoj060603/recomb2027-intervention-codes
