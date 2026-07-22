# EntangledUnits refactor — plan (scratch, not for commit)

Branch: `entangled-refactor`. Scratch working notes; delete before merging.

## Goal

Drastically simplify handling of "entangled units". End state: **`EntangledSystem`
disappears; its capability is absorbed into a generalized `System` that may be
entangled (`nparts > 1`).** The `EntangledUnits/` wrappers collapse.

## Status

Done:
- `MeasureSpec` carries an `nparts` axis.
- `SpinWaveTheory(esys)` → thin wrapper
- `SampledCorrelations(esys)` → thin wrapper

**Net: the calculators (SWT, SC) are ready. What remains is the `System` side.**

- Phase A done: design pinned (see Phase A section). `Entanglement` struct wrapping a
  full `bare_system::System` + mapping/dynamics metadata; `entanglement` optional field
  on `System`. Eager sync + full bare_system kept for B; slimming deferred to B′.

### Transitional seam introduced this session
`SampledCorrelations(sys; measure, ..., crystal, origin_crystal)` — the entangled path
passes `crystal = esys.sys_origin.crystal` (physical crystal to report against, while
coherents are sampled from `esys.sys`). This explicit `crystal` kwarg is the one residual
coupling. It gets removed in Phase D once `System` carries physical geometry itself.

## Key structural facts (why this is tractable)

- `EntangledSystem` (`src/EntangledUnits/TypesAndAliasing.jl:37`) holds:
  - `sys` — contracted SU(N) System, `nunits` sites/cell (one coherent per unit).
  - `sys_origin` — physical System, `natoms_origin` sites/cell (real dipoles).
  - `contraction_info :: CrystalContractionInfo` — `forward[atom]=(unit,k)`,
    `inverse[unit]=[InverseData(site, offset)]`. `offset = physical_pos - unit_pos`.
  - `source_idcs`, `dipole_operators` — dynamics metadata / cached product-space spin ops.
- Correlations stay **flat** over physical positions (`npos = nparts·nunits`);
  the `nparts` structure is bookkeeping (offsets, form factors, which coherent to read),
  not a correlation-tensor axis. Confirmed: `positions[p] = unit_pos + offset[k,unit]`.
- Retrieval only needs the physical *crystal* + positions, NOT a full second System.
  → biggest simplification lever: `sys_origin` may collapse to "physical crystal +
  dipoles buffer" rather than a full `System`.

## Phases (dependency-ordered; each independently landable + golden-gated)

### Phase A — Characterize residual coupling (research) — DONE
Inventory what `EntangledSystem` provides beyond a plain `System`:
1. Two-system pairing (contracted `sys` ↔ physical `sys_origin`) + `contraction_info`.
2. Dynamics sync: `step!`/`randomize_spins!`/`minimize_energy!` act on `sys`, then
   `set_expected_dipoles_of_entangled_system!` writes physical dipoles into `sys_origin`
   (`src/EntangledUnits/EntangledUnits.jl:342`).
3. Construction: `entangle_system` builds contracted interactions from physical bonds
   (`src/EntangledUnits/EntangledUnits.jl:243`).
4. Aliasing: ~20 pass-through methods in `TypesAndAliasing.jl`.

**Conclusions:**
- `sys_origin` (physical System) is *read* only for: `crystal` (+`root`/`origin` for
  reshaping), `Ns`, `gs`, and the `dipoles` buffer. Never for `interactions_union`,
  `coherents`, `params`, `ewald`, `κs`, `rng` — those live on the contracted `sys`.
- So it **could** collapse to a slim `{crystal, Ns, gs, dipoles}` record. **But** four
  capabilities call `System`-level methods on it (reshaping trio + `position_at`/
  `position_to_site`; `ssf_custom(f, ·)`; `set_field!`; `magnetic_moments`/`plot_spins`).
  Collapsing to a slim record forces crystal-level rewrites of those. **Decision: keep a
  full `System` (renamed `bare_system`) in Phase B; slimming is a later, independent
  Phase B′ — do NOT bundle it into the riskiest step.**
- Physical `dipoles` buffer is write-only w.r.t. SC/SWT retrieval (they read
  `sys.coherents` via product-space ops); only `magnetic_moments`/`plot_spins` read it.
  Lazy sync is possible (kills the "meaningless sync" TODO at EntangledUnits.jl:304) but
  **deferred** — Phase B keeps eager sync verbatim to minimize risk.

**Pinned design (dedicated struct, `entanglement::Union{Nothing,Entanglement}` on System):**
```julia
struct Entanglement
    bare_system      :: System                                 # physical (uncontracted) system; Phase B′ may slim
    contraction_info :: CrystalContractionInfo
    source_idcs      :: Array{Int64, 4}
    dipole_operators :: Vector{NTuple{3, Matrix{ComplexF64}}}
end
```
`is_entangled(sys) = !isnothing(sys.entanglement)`. The unified `System` *is* today's
`esys.sys` (contracted) + `.entanglement` back-pointer. `EntangledSystem(sys, units)`
becomes a constructor returning such a `System`.

Method → unified home:
| Today (`esys.…`)                                   | Unified `System`                                  |
|----------------------------------------------------|---------------------------------------------------|
| `esys.sys.*` (coherents, interactions, buffers)    | `sys.*` directly                                  |
| `esys.sys_origin`                                  | `sys.entanglement.bare_system`                    |
| `contraction_info`/`source_idcs`/`dipole_operators`| `sys.entanglement.*`                              |
| `eachsite(esys)`→`eachsite(sys_origin)`            | `eachsite` dispatch: entangled → iterate bare_system |
| `eachunit(esys)`→`eachsite(sys)`                   | `eachsite(sys)` (units are native sites)          |
| `orig_crystal`,`energy_per_site`,`magnetic_moments`,`plot_spins`,`set_field!`,`set_coherent!` | same bodies, read `sys`+`entanglement.bare_system` |
| dynamics sync                                      | Phase C: fold into `step!(sys)` when `is_entangled`|
| `ssf_custom(f,esys)`→`entangled_measure`           | `ssf_custom(f,sys)` branches on `is_entangled`    |
| reshaping trio                                     | operate on `sys`+`bare_system`, rebuild `contraction_info` |

Phase B plumbing touch-points: `clone_system` (clone inner `bare_system`), positional
`System(...)` constructor calls at System.jl:88 & 161 (pass `nothing`), `show`.
Golden invariant: "Dimer Tests" intensities stay bit-for-bit (B only relocates state).

### Phase B — Give `System` an optional entanglement mapping (PIVOTAL, riskiest)
Add `entanglement::Union{Nothing, Entanglement}` as an **optional field on `System`**
(`nothing` for ordinary). Ordinary = `nothing` path. An entangled `System` = today's
`EntangledSystem.sys` + `.entanglement` back-pointer holding a full `bare_system::System`.
Struct shape and method→home mapping pinned in the Phase A section above. Scope B to pure
state relocation: keep eager dipole sync and the full `bare_system` (no slimming — that's
B′). Golden-gate on "Dimer Tests".

**Phase B landed:**
- `abstract type AbstractEntanglement` declared in `System/Types.jl` (breaks recursive
  `System`↔`Entanglement` dependency). New non-const field `entanglement::Union{Nothing,
  AbstractEntanglement}` appended to `System{N}`.
- Concrete `struct Entanglement <: AbstractEntanglement` in `EntangledUnits/
  TypesAndAliasing.jl`: `{bare_system::System, contraction_info, source_idcs,
  dipole_operators}`.
- `is_entangled(sys) = !isnothing(sys.entanglement)`; `clone_entanglement` (nothing
  fallback in System.jl, Entanglement method in TypesAndAliasing.jl).
- Field populated in the `EntangledSystem(sys, sys_origin, contraction_info)` funnel; the
  invariant `sys.entanglement.bare_system === esys.sys_origin` is re-established in
  `clone_system(::EntangledSystem)` too.
- 3 positional `System(...)` call sites (System.jl:88,161; Reshaping.jl:118) pass
  `nothing`. Nothing yet *reads* `sys.entanglement` — that's Phase C.
- Verified: "Dimer Tests" 14/14 bit-for-bit; contraction+reshaping group 36/36; live
  invariant checks (ordinary→nothing/false; entangled→populated & aliased; clone deep +
  invariant preserved).

### Phase C — Fold dynamics + construction onto unified `System`
- `step!(sys)` on entangled `System` does coherent evolution + (lazy) dipole expectation
  internally; the sync becomes an internal detail, not a wrapper obligation.
- `EntangledSystem(sys, units)` → constructor/`entangle!` returning a `System` with the
  mapping field populated.

### Phase D — Delete `EntangledSystem` + wrappers
Wrappers (`EntangledSampledCorrelations`, `EntangledSpinWaveTheory`) and most of
`TypesAndAliasing.jl` evaporate. **Remove the transitional `crystal` kwarg from
`SampledCorrelations`** — it reads physical geometry from the system's mapping field.

Ordering rationale: B is the linchpin; A de-risks it cheaply; C/D mechanical once B lands.
Do NOT attempt B+C+D in one session. Golden invariant throughout: "Dimer Tests" golden
intensities in `test/test_entangled_units.jl`.

## Fresh-session kickoff prompt (Phase A)

> Continue the EntangledUnits refactor. SC and SWT are already thin wrappers consuming
> `nparts` MeasureSpecs. Next: Phase A — inventory exactly what `EntangledSystem`
> provides beyond a plain `System` (two-system pairing, dynamics sync, construction,
> aliasing) and map each to its future home on a unified `System`. Read
> `src/EntangledUnits/TypesAndAliasing.jl` and `EntangledUnits.jl`. No code changes yet
> — produce the inventory and confirm whether `sys_origin` can collapse to
> crystal+buffer. Test recipe is in CLAUDE.local.md.

## Notes
- Start Phase B+ in a FRESH julia-mcp session (struct redefs; Revise can't track).
- Commit the current SC work before starting the next phase (clean `git status`).
- Test recipe: see "Tests" section of CLAUDE.local.md (run_tests at repo ROOT via
  julia-mcp; LOAD_PATH trick; don't Pkg.develop into test env).
