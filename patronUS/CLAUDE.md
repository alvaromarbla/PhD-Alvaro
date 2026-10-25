# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

PATRONUS (Performance Analysis, Trajectory Routing and Optimized Navigation for Unmanned Systems) is a MATLAB tool for performance optimization of a tilt-rotor VTOL UAV, developed as part of a PhD at the University of Seville. The git root is the parent `scripting/` folder; active code lives in `patronUS/`. `../legacy/` holds the pre-refactor scripts (fmincon-based, plus a vendored OptimTraj) and is kept only for reference.

## Running

There is no build, lint or test suite. Everything runs inside MATLAB:

- Requires CasADi for MATLAB on the path (`import casadi.*`, IPOPT solver). Nothing in the repo calls `addpath`, so `patronUS/` and all its subfolders must be added to the MATLAB path manually (e.g. `addpath(genpath('patronUS'))`).
- Entry point: run `main_solvers/main_cruise_optimizer.m` (a script, not a function). The `versions` struct at its top selects aircraft, model variants, constraint set and objective.
- `main_solvers/main_cruise.m` is an older fmincon/SQP version that no longer matches the current `cruise_config` signature.

## Architecture

The optimizer is a single-point (steady cruise) NLP built with CasADi `Opti` and solved with IPOPT. The flow in `main_cruise_optimizer.m`:

1. **`versions` struct** chooses `AC_version` (`"1Wing"` / `"2Wings"`), per-coefficient model variants (`CL/CD/CM/CT/CP/CH_version`: `'analytical'`, `'dummy'`, `'bspline'`, `'zeros'` depending on the model), `constraints` and `objective`.
2. **`config/cruise_config.m`** returns `params` (mass, geometry, propeller data, atmosphere), `bounds.vars` and `models`. Each field of `bounds.vars` is `struct('min','max','init')` and **defines a decision variable**: the main script builds the `X` vector dynamically from `fieldnames(bounds.vars)`, so adding/removing/commenting a variable there changes the problem. Variables are accessed by name via `X_struct.<name>` in objectives/constraints.
3. **`main_solvers/models_config.m`** builds `models`, a struct of CasADi `Function`/`interpolant` objects (`models.CL_wing_lookup(alpha)`, `models.CT_lookup(J, phi)`, …) from the factories in `vehicle/aero/` and `vehicle/prop/`. Each factory switches on its version string. Propeller coefficients take advance ratio `J = V/(n·D)` and inflow angle `phi = alpha + epsilon`. Models must stay CasADi-symbolic (use smooth approximations like `sqrt(x^2+1e-6)` instead of `abs`).
4. **`config/patronus_registry.m`** maps the `versions.constraints` / `versions.objective` keys to function handles in `segment/cruise/`. New formulations are registered there, not in the main script. `validate_config.m` checks the keys exist.
   - Constraint functions: `[c, ceq, dbg] = f(X_struct, models, params, bounds)`, with `c <= 0` and `ceq == 0`. `dbg` is an optional struct of intermediate expressions used for post-processing.
   - Objective functions: `obj = f(X_struct, models, params)`, minimized and non-dimensionalized.
5. **Post-processing**: for `2Wings`, `segment/cruise/cruise_postp_2Wings_dbg.m` evaluates the `dbg` struct at the solution and `generate_results_page_2Wings.m` plots it. For `1Wing`, `main_solvers/cruise_postp_1Wing.m` → `utilities/generate_results_page_1Wing.m`.
6. **IPOPT diagnostics**: `IterHistory` (handle class) records `X` and `g` through `opti.callback`. Then `computeConstraintViolation`, `computeBoundSlack`, and `utilities/plotIpoptconvergence`, `plotVariableEvolution` and `plotConvergenceMetrics` plot convergence. On a failed solve the script falls back to `opti.debug.value` and skips the results page, but still runs the diagnostics.

### Aircraft variants

- **`1Wing`**: 2-DOF point-mass force balance (`cruise_constraints_2GDL`, or `cruise_constraints_nacelle` for `1Wing_Nacelle`). Variables are `V, gamma, alpha, epsilon, n`, with a single CL/CD model for the whole aircraft.
- **`2Wings`**: 3-DOF longitudinal trim (`cruise_constraints_3GDL`). Wing, fuselage and tail aero are modelled separately, with tail downwash from `downwash_calc`. It has two tilting propeller groups (`epsilon1/2`, `n1/2`), elevator `deltae` and a nacelle drag model that solves the induced-velocity quartic by fixed Newton iterations. It also includes pitching-moment balance using `params.geo` arm lengths. Variables use `theta` (with `alpha = theta - gamma`) rather than `alpha`.

## Gotchas

- The main script always calls the constraint function with three outputs (`[c, ceq, dbg]`). The 1-Wing constraint functions currently return only `[c, ceq]`, so they need a `dbg` output to run through the current main script.
- `utilities/plotforcesangles_imag.m` loads `fuselage_sideview.png` and `nacelle_sideview.png`, which are not in the repo; they must be on the MATLAB path.
- Figure exports write vector PDFs to `patronUS/Images/` (gitignored), resolved relative to the source file via `mfilename('fullpath')`. Keep paths machine-independent: the project is worked on from several machines.
- Comments and some names mix English and Spanish, and "GDL" means DOF (grados de libertad).
