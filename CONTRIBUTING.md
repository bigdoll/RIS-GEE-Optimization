# Contributing

Thanks for your interest in improving this repository. It accompanies an IEEE
Transactions on Communications paper, so the priority is that the code stays a
faithful, readable implementation of the published algorithms.

## Ground rules

- **Keep the math traceable.** Every function maps to an equation or algorithm
  in the paper (see [`docs/equation_mapping.md`](docs/equation_mapping.md)). If
  you add or change a function, update its header comment and that table.
- **Don't silently change numerics.** Function names follow the paper's
  notation and the computational lines are intentionally preserved. If a change
  alters results, say so in the PR description and explain why.
- **One concern per pull request.** Small, focused PRs are easier to review.

## Project layout

```
src/common/        shared system model (channels, noise, power, func_R)
src/active_ris/    active RIS    -> approach1_alternating / approach2_mmse
src/passive_ris/   passive RIS   -> approach1_alternating / approach2_mmse
examples/          runnable entry points
tests/             integrity checks
docs/              system model & equation mapping
```

Within each approach, subproblems live in `ris_reflection/` (Algorithm 1 / 4)
and `power/` (Algorithm 2 / 5); the driver (Algorithm 3 / 6) sits one level up.

## Before opening a PR

1. **Static check (no MATLAB needed):**
   ```bash
   python3 tests/check_dependencies.py
   ```
   This must pass: it confirms there are no orphan functions and that every
   intra-repository call resolves.

2. **MATLAB smoke test (if you have MATLAB):**
   ```bash
   matlab -batch "addpath tests; smoke_check"
   ```
   Runs without a solver and exercises the shared pipeline on a tiny instance.

3. **Full run (if you have MATLAB + CVX + MOSEK):** run
   `examples/run_active_ris.m` and `examples/run_passive_ris.m` and confirm the
   curves still look sensible.

## Coding style

- MATLAB function help block at the top of every file (purpose + paper ref).
- Prefer descriptive comments over renaming paper-notation variables.
- Keep CVX models inside `*cvxopt*` files; the solver is selected with
  `cvx_solver mosek` (swap to `sedumi` if you don't have a MOSEK license).

## Reporting issues

Open a GitHub issue describing the scenario, parameters used, MATLAB/CVX/MOSEK
versions, and the observed vs. expected behavior.
