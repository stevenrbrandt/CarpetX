# ODESolvers

| Author(s)      | Erik Schnetter and Liwei Ji |
|:---------------|:----------------------------|
| Maintainer(s)  | Erik Schnetter and Liwei Ji |
| Licence        | LGPL |


## Purpose

Solve systems of coupled ordinary differential equations


## Tag requirements

A group with an `rhs=` tag is integrated by ODESolvers. ODESolvers checks its tags at startup (`CCTK_PARAMCHECK`) and aborts with an error naming the group if they are not allowed:

- It must be checkpointed, i.e. it must not set `checkpoint="no"`, whatever its `evolve` tag. A checkpoint would omit its state and, under subcycling, its Runge-Kutta source bands, so recovery could not restore it. Enforced by `ODESolvers_CheckTags`, which tests this rule first.
- It must be evolved, i.e. it must not set `evolve="no"`. Regrid would not reliably refill its state. Enforced by `ODESolvers_CheckTags`.
- Under subcycling (`CarpetX::use_subcycling = yes`) it must also have a single timelevel. Enforced by `ODESolvers_CheckTimelevels`.

The default tags (no `checkpoint` or `evolve` tag) satisfy the first two rules.

## RK4-2 and RK4-3

`RK4-2` and `RK4-3` are the multistep Runge-Kutta methods from
arXiv:2603.05763. They keep previous RHS evaluations in time levels of the
RHS groups. `ODESolvers::rhs_timelevels = 0` allocates the minimum the
current method needs: 2 for `RK4-2`, `RK42-IMEX`, and `RK42-IMEX2L`, 3 for
`RK4-3`, and 1 for every one-step method. Each RHS group still declares
`TIMELEVELS=4` in its `interface.ccl`, which is the ceiling, because
`ODESolvers::method` is steerable. A positive `rhs_timelevels` must be at
least that minimum and at most the declaration.
`RK4_dash_2_sol` selects published solution (1) (the default) or (2). The
first steps, and any step whose history was invalidated by regridding or
recovery, take a classic RK4 step to fill that history. These methods are
not available when `CarpetX::use_subcycling` is `yes`.

`RK42-IMEX` is the same RK4-2(1) tableau on the explicit RHS (fluxes, and
whatever else `ODESolvers_RHS` evaluates), with one backward-Euler source
step after each positive stage abscissa and after the final update. The
source step is `ODESolvers_ImplicitStep`, which nuX implements in
`nuX_M1_CalcUpdate`. Solution (2) is rejected: its first abscissa is
negative. A negative abscissa keeps the explicit update and, when
`nuX_Base::rF` is an evolved group, puts the momentum back to the last
relaxed stage. The step that fills a missing history is IMEX Euler, not
classic RK4, because the explicit RHS does not contain the collision term.
This is not an order-4 implicit method. Its source step is backward Euler,
so the source part is first order.

`RK42-IMEX2L` keeps those same RK4-2(1) flux weights and replaces the
source step with a two-stage L-stable DIRK (`gamma = 1 - sqrt(2)/2`).
Both implicit solves use step `gamma*dt`. The source tendency from the
first solve is added into the later flux stages, including the negative
abscissa, which is a weight and not a solve. There is no momentum
restore. The finished step is the second solve. A missing history takes
one classic RK4 flux step coupled to that same DIRK. Later steps are
RK4-2 on the flux. Those steps match RK4 and RK4-2 when the implicit
kernel leaves its saved base unchanged, so the recovered source is zero.
nuX floors and projects the radiation in the pre-implicit PostStep for
this method only, which keeps that repair out of the source tendency.
IMEX42L instead rebuilds the step from the recovered source, so the same
repair stays inside its source. The source part is order 2. This is not
an order-4 implicit method, and subcycling is refused.

## IMEX32L

`IMEX32L` is the IMEX-SSP3(3,3,2) scheme in Table 5 of Pareschi and
Russo, "Implicit-Explicit Runge-Kutta Schemes and Applications to
Hyperbolic Systems with Relaxation," Journal of Scientific Computing
25 (2005) 129–155
(<https://doi.org/10.1007/s10915-004-4636-4>, arXiv:1009.2757). The
explicit part is SSPRK3, with nodes \(0\), \(1\), and \(1/2\). The
implicit part is an L-stable DIRK with diagonal
\(\gamma = 1 - 1/\sqrt{2}\) on every stage, so each step takes three
implicit solves of size \(\gamma\Delta t\). The implicit nodes are
\(\gamma\), \(1-\gamma\), and \(1/2\). The flux and the source of one
stage are evaluated at those two different times. The weights are
\(1/6\), \(1/6\), and \(2/3\). The finished step is that weighted
combination, so a repair made inside the implicit solve stays in the
recovered source. The combined scheme is order 2. With the source absent
it is SSPRK3, which is order 3.

## IMEX42L

`IMEX42L` is the IMEX-SSP3(4,3,3) scheme in Table 6 of Pareschi and
Russo, the same paper as `IMEX32L`. The explicit part is SSPRK3, with
a leading stage whose flux weight is zero. Its nodes are \(0\), \(0\),
\(1\), and \(1/2\). The implicit part is an L-stable DIRK. Every
diagonal is the Table 6 value \(\alpha = 0.24169426078821\), so each
step takes four implicit solves of size \(\alpha\Delta t\). The
implicit nodes are \(\alpha\), \(0\), \(1\), and \(1/2\). The weights
are \(0\), \(1/6\), \(1/6\), and \(2/3\). The first stage solves
\(Y = y_n + \alpha\Delta t\, g(Y)\). It does not evaluate the stiff
source explicitly.

The scheme that used to be stored under this name paired classic RK4
with an implicit matrix whose first diagonal was 0. That matrix is the
I42L(4,2) tableau. Its first stage stored the stiff rate, and the next
stage added \(\Delta t/4\) times that rate before primitive recovery.
Table 6 replaces it. The keyword is unchanged. The flux weights in the
Runge-Kutta method table are the Table 6 weights, not the RK4 weights.

The second implicit stage begins with \(-\alpha\). A later stage
predictor can leave the range of the data. The predictor is a
combination of solved stages. The combined scheme is order 3. With the
source absent it is SSPRK3, which is order 3. A repair made inside an
implicit solve stays in the recovered source, as it does for `IMEX32L`.

## Subcycling

Add the following parameters to your parameter file

```
CarpetX::use_subcycling = yes
CarpetX::restrict_during_sync = no
```

Memory: under subcycling each refined level holds, per evolved group, two extra zero-ghost work buffers over its coarse-fine boundary footprint, one at the coarse and one at the fine resolution (the latter has about 8 times the cells of the former). They come on top of the `num_rk_stages + 1` coarse-resolution source bands that each refined level keeps per evolved group, which hold its parent's start-of-step state and stage derivatives on the same footprint. The bands and the two work buffers are allocated together, when the parent level first stores its old state after the refined level was made, and stay allocated until the next regrid remakes the level, so that the per-stage boundary fill allocates nothing. The bands of a group cover the footprint of that group's own prolongation stencil, so groups with different `prolongation_type` / `prolongation_order` tags may be mixed freely. Runs without subcycling never allocate any of these.


## IMEX methods

Both drivers use `src/imex_tableaux.hxx`. Implicit increments are recovered
before `PostStep` numerical projections. No extra final implicit projection
is applied. Flux registers use the same explicit weights as the accepted state.

- `IMEX42L` is Table 6 of Pareschi and Russo, described above. It replaces
  Table 10 of [Izquierdo et al.](https://arxiv.org/html/2211.00027v2). The
  explicit part is SSPRK3 and the additive scheme is order 3. Its implicit
  part is L-stable. That alone is not a stiff coupled-transport guarantee.
- `IMEX32L` now uses the L-stable IMEX-SSP3(3,3,2) in Table 5 of
  [Pareschi and Russo](https://arxiv.org/pdf/1009.2757). Its explicit part is
  SSPRK3, and its implicit diagonal is `gamma=1-1/sqrt(2)` at all three stages.
  Explicit times `(0,1,1/2)` differ from implicit times `(gamma,1-gamma,1/2)`;
  each partition's time-dependent fields and boundaries are refreshed there.
  **This deliberately changes results from the former IMEX32L.** The old
  printed CK tableau has `b.A_I.1=5/6`, not `1/2`, and is only first order
  for pure implicit decay. Matching the printed table was not sufficient.
- `IMEX122` is explicit/implicit midpoint. The second RHS is evaluated at the
  midpoint and the implicit increment is scaled exactly once. It is second
  order and A-stable, not L-stable, and is not supported under subcycling.
  `Implicit Euler` remains a deprecated compatibility alias with a startup
  warning: it does not mean backward Euler.

Under subcycling, IMEX coarse-fine data use accepted parent endpoint secants.
Both partition times, including the first nonzero implicit stage, are evaluated
at their actual parent-step fractions. This is linear time interpolation, not
a new high-order IMEX dense-output formula.

For the scalar split `y'=(lambda_E+lambda_I)*y`, at fixed `z_E=lambda_E*dt`,
the unprepared stiff limit `z_I -> -infinity` is zero for corrected `IMEX32L`
and for Table 6 `IMEX42L`. The replaced Izquierdo Table 10 limit was `-z_E/2`.
Midpoint's pure implicit amplification tends to
`-1`. Do not equate pure-implicit L-stability, nonautonomous order, prepared
diffusion accuracy, or nonlinear radiation-matter stability: test each.

Corrected `IMEX32L` is not globally stiffly accurate. In a relaxation problem,
its internal stages approach the equilibrium manifold and its conserved
variables approach the explicit reduced-system update, but the accepted fast
variables need not lie exactly on that manifold. The next first implicit stage
relaxes them again. This distinction matters when examining output flux
speckles; it is not permission to append an unweighted implicit projection.

Run the source-reading math checks without a Cactus build:

```
python3.11 ODESolvers/test/test_imex_tableaux.py
```

They check the actual shared coefficients, additive order, scalar autonomous
and time-dependent convergence, stiff limits, and partition-time/source
contracts; all 12 currently pass. They do not execute the compiled schedules.
The paired nuX/ODESolvers CUDA build, short source-only relaxation runs and
a four-rank/four-GPU two-step IMEX42L TOV smoke have also completed. These
checks do not establish all-method AMR order, nonlinear diffusion accuracy
or long coupled TOV stability; those validations remain required.
