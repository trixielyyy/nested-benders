# Nested Benders for the three-stage multiple knapsack

Status: **prototype, not executed.** The code builds. No run is verified. The
extensive sweep holds the Gurobi WLS license, so this directory starts no run.

## 1. Can nested Benders solve this problem?

The question is about convexity, not about the code.

Nested Benders works when each node value function is convex in the parent
decisions. Cuts come from duals. Duals describe a convex function.

In this benchmark the first two stages are binary. So:

- The root master holds `x, w` (binary) and `theta_lambda`.
- Node `lambda` holds `y^lambda` (binary) and `psi_s`.
- The stage-3 subproblem is an LP in `z^s`.

The node value `V_lambda(x)` is the minimum of a mixed-integer problem. That
function is **not convex** in `x`, so it has no dual and no supporting
half-space. A Benders cut in the usual sense does not exist.

A cut built from the node's LP relaxation is a different object. That cut is
valid: `V^LP` is convex, the hyperplane lies below `V^LP` everywhere, and
`V >= V^LP`. Validity holds. Certification does not, because the cut prices the
node at its relaxation value. The master can then choose a point whose true
integer cost is higher than its estimate.

This is a known limit of Benders applied to integer subproblems. It is not a
defect of the implementation.

Three options are open.

**Option A: relax stage-2 integrality.** The node becomes an LP. The value
function becomes convex and the cuts become valid. The bound converges to the
LP relaxation, not to the MIP optimum. The prototype implements this option.

**Option B: keep `y` in the master.** This is standard Benders, and
`knapsack/` already runs it. The nesting disappears, because no integer
information reaches a node.

**Option C: integer cuts.** Option C is exact, with an exponential worst case.
Section 2 covers the variants.

## 2. How to choose the integer cuts

Three schemes appear in the literature. They differ in what they remove and
what they cost.

**Integer optimality cut (Laporte and Louveaux).**
`theta >= V* - M (1 - delta)`, where `delta` is a binary that equals one only
at the incumbent pair `(x, y)`. The cut is valid and exact. It needs one
`delta` binary and one big-M value per cut. It converges slowly. The big-M
value carries numerical risk when objective ranges are large.

**No-good cut on `y`.**
The cut forbids the current node pattern:
`sum_{j in y+} y_j - sum_{j in y-} y_j <= |y+| - 1`.
It is valid, exact, and needs no big-M. It also carries no cost information,
so the root master learns nothing from it. The outer bound stalls while the
node enumerates patterns. In the worst case the node visits `2^|y|` patterns.

**Branch and Benders at the node.**
The stage-2 variables enter a branch-and-bound tree. Stage-3 cuts are added
only at integer-feasible node solutions, as lazy constraints. Each cut stays
valid, because `y` is fixed when that cut is generated. The existing
`ScenarioBenders` class already works this way inside one node.

**Recommendation.** Use the branch-and-Benders variant for exactness, or keep
option B and skip nested Benders. Use option A only to produce a lower bound.

No scheme above rescues the root cut. A cut over `x` is valid only if the node
is solved for all `x`. Any exact scheme therefore also needs a branch-and-bound
over `x`. That is integer L-shaped, not nested Benders.

## 3. Code reuse

**Reused as-is:**

- `knapsack/ProblemParams.cpp` parses `.3stage.cor`, `.3stage.tim` and
  `_extended.sto`. That is about 250 lines nobody needs to rewrite.
- `knapsack/GurobiSolve3Stage.cpp` builds the full extensive form. The
  prototype calls it for the reference optimum.

**Not reused:**

- `BendersMaster3Stage` and `ScenarioBenders` assume one two-level master. The
  node structure is different. The record `ScenarioOptimalityCut` has the same
  fields as `Stage3Cut` in the prototype, so the two formats match by design.

Build rule: the prototype includes `GurobiSolve3Stage.h` and links its `.cpp`.
No file in `knapsack/` changes.

## 4. What the prototype does

`nested_benders.cpp` builds three model families:

1. Root master: `x, w` binary, `theta_lambda` free, stage-1 rows.
2. One node per `lambda`: `y^lambda`, `psi_s` free, stage-2 rows with `x` moved into
   the right-hand side.
3. One subproblem per `s`: stage-3 LP. The dual point comes from `Pi`.

Two modes select the node type and the root cut.

**Exact mode (default).** The node is a MIP. The inner loop solves the node, checks
every stage-3 subproblem at the node solution, adds the violated cuts, and solves the
node again. That is Benders on the node, with the node master carrying integer `y`.
When no cut is violated, the node value `V_lambda(x_hat)` is exact.

The root cut is an integer optimality cut. Each cut adds one binary `delta`, forces
`delta = 1` only when `x = x_hat`, and states `theta_lambda >= V_lambda(x_hat)` under
that indicator. The cut claims nothing at other `x`, so it needs no convexity.

**Relaxed mode (`--relaxed`).** The node is an LP. The root cut is the usual
subgradient cut. The bound converges to the LP relaxation.

### Why the exact mode terminates
The root value is a lower bound on the optimum. That needs one property: `theta_lambda`
has a finite lower bound. Here the bound is zero, which holds because every cost
coefficient is nonnegative. Without it the root gives every unvisited pattern `-inf`
and no test can fire.

The cut is the Laporte and Louveaux (1993) integer optimality cut in its big-M form:

```
theta_lambda >= V_lambda(x_hat) - M * d(x, x_hat)
```

where `d(x, x_hat)` counts the stage-1 coordinates that differ. With binary `x`, that
count is linear: `sum_{x_hat_j=1} (1 - x_j) + sum_{x_hat_j=0} x_j`. The cut needs no
auxiliary binary.

The value `M = V_lambda(x_hat)` is valid here. Every cost is nonnegative, so
`V >= 0` and `theta_lambda >= 0`. One deviation then makes the row read
`theta_lambda >= 0`, which is already the variable bound. The run refuses exact mode
when a cost coefficient is negative, because `M = V` stops being a valid bound there.

The loop stops when the root lower bound reaches the true cost of the current pattern.
Each iteration adds one cut, so the loop terminates. The worst case visits every
stage-1 pattern, which is exponential in the number of stage-1 variables.

### Why not a no-good cut
A no-good cut needs no new binary either:

```
sum_{x_hat_j=1} x_j - sum_{x_hat_j=0} x_j <= |x_hat| - 1
```

It removes the current pattern and says nothing about `theta`. The root then pays a
cost of zero for every unvisited pattern, so it learns no cost information. The method
degrades to enumeration of stage-1 patterns. The integer optimality cut carries the
cost `V` at the visited pattern, so the bound rises with each iteration and the search
can prune. The no-good cut is dominated here, and the code does not implement it.

### What the run prints

The reference MIP and LP values, the nested objective, the certification flag, the
number of node and subproblem solves, the size of the root master, the size of one
node master, and both gaps. The size lines give the memory argument for the nested
structure: the root carries no stage-2 or stage-3 variables.

## 5. Verification plan (not done)

1. Check the extensive builder. Its MIP value must equal the value in
   `results_rvar/B28/<config>_smkp_<n>.extensive.log`.
2. Run the exact mode. The reported objective must equal the extensive MIP value,
   and the run must print `(certified)`.
3. Run the relaxed mode. The reported bound must stay at or below the MIP optimum
   and at or above the LP relaxation.
4. Record outer iterations, node solves, and master sizes for each configuration.

Steps 1 to 3 need a license slot. Run them when the sweep is idle.

## 6. Known gaps

- Feasibility cuts are absent. An infeasible node or subproblem aborts the run.
- The inner node loop stops after `maxInner = 40` solves. A cap hit sets a flag, and
  the run reports that the node value is not confirmed.
- The exact mode adds `n1 + 2` rows and no binary per node per outer iteration.
  Memory grows with the iteration count.
- No warm start applies between outer iterations.
- The relaxed subgradient assumes Gurobi's sign convention on `Pi`. Step 3 tests that
  assumption.
