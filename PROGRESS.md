# Progress log

State: the prototype builds and is **not run**. The diagram is generated,
checked and ready to share. The repository holds both.

## 2026-11-03 — repository

- Created `github.com/trixielyyy/nested-benders`. The repository is private.
- First commit `7b4d244` on `main` holds six files: `README.md`, `Makefile`,
  `nested_benders.cpp`, `nested-benders-flow.html`,
  `.archify/workflow-nested-benders-20261103-200233/candidate.json` and
  `.gitignore`.
- Commit identity: `trixielyyy <trixielyyy@users.noreply.github.com>`. It
  matches the manuscript repository.
- `.gitignore` keeps the compiled binary and the archify receipts out of git.

## 2026-11-03 — prototype

`nested_benders.cpp`, 502 lines.

- Three model families: root master over stage-1 patterns, one integer node
  master per stage-2 scenario, one stage-3 LP per third-stage scenario.
- Exact mode is the default. The node runs as a MIP. Stage-3 cuts enter the
  node between solves, so the node value is exact at the tested pattern.
- The root cut is the Laporte and Louveaux integer optimality cut in big-M
  form, `theta >= V - M d(x, x_hat)`, with `M = V`.
- `theta` starts at 0. This bound is valid because every cost coefficient is
  nonnegative, and the run stops when the root bound meets the true cost of the
  current pattern.
- `--relaxed` keeps the old lower-bound experiment.
- A stage-1 variable that is not binary forces `--relaxed`. The code reports
  it.

Three defects were found and fixed during the build:

1. The first stop test compared the root value against the node sum. That test
   is invalid when `theta` is free at unvisited patterns.
2. The `delta` link rows were one-sided, so the root could ignore the cut. The
   `delta` construction was then dropped for the big-M form.
3. `theta` had no lower bound, so no test could fire.

Build: `make` completes with no errors and no warnings.

## 2026-11-03 — diagram

Generated with archify 3.2.0, diagram type `workflow`.

- Legend fix. Archify filled the legend from its own example, so it read
  "User UI", "Agent logic", "Policy", "Tool action". It now reads Instance,
  Solver step, Benders cut, Bound check.
- unslop pass, five changes:
  - preset `blueprint` dark to `editorial`.
  - colour collapse from six bright strokes to three roles.
  - motion off.
  - title and subtitle tightened, lane labels shortened.
  - three card dots reduced to one slate colour.
- Gates `validate`, `deliver`, `check` and `browser-check` pass.
- Artifact: 846,692 bytes, sha256 `7e71cf0…`.

## Open items

1. Run the prototype once. This needs a Gurobi licence slot while the knapsack
   sweep is idle.
2. Check the extensive builder against
   `knapsack/results_rvar/B28/<config>_smkp_<n>.extensive.log`.
3. Check the nested bound against the LP relaxation.
4. Decide how to host the page. The repository is private, so GitHub Pages
   needs a public repository or a paid plan.
5. Remove the remaining viewer chrome if wanted: the `01 /` lane numbering,
   the " Diagram" title suffix, the pulse dot and the toolbar. All four live in
   the archify template. A hand patch ends regeneration.
6. Commit and push the HTML and the candidate. Both are modified in the working
   tree and neither is pushed.

## Context outside this directory

The knapsack extensive sweep stood at 162 of 180 runs with a result at
2026-11-03 19:45. Nothing in this directory touched it.
