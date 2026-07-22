# Reshaping of an entangled `System`. Both the contracted `sys` and its physical
# `bare_system` are reshaped by the ordinary (non-entangled) backbone, then the
# entanglement mapping is rebuilt geometrically from the immutable `units` truth
# (see `rebuild_entanglement!`). Because the rebuild starts from the single-cell
# truth rather than transforming the previous mapping, units may straddle cell
# boundaries of the reshaped system. These are dispatched from the ordinary
# `reshape_supercell`/`repeat_periodically`/`resize_supercell` (Reshaping.jl).

# Apply the ordinary reshaping `op` (a closure of one system argument) to both
# the contracted `sys` and its physical `bare_system`, then rebuild the mapping.
function reshape_entangled(sys::System, op)
    (; bare_system, units) = get_entanglement(sys)
    sys_new = reshape_contracted(sys, op)
    bare_new = op(bare_system)
    return rebuild_entanglement!(sys_new, bare_new, units)
end

reshape_supercell_entangled(sys::System, shape) =
    reshape_entangled(sys, s -> reshape_supercell(s, shape))

repeat_periodically_entangled(sys::System, counts) =
    reshape_entangled(sys, s -> repeat_periodically(s, counts))

resize_supercell_entangled(sys::System, dims::NTuple{3,Int}) =
    reshape_supercell_entangled(sys, diagm(Vec3(dims)))

# Flatten an entangled `sys` into a single `(1,1,1)` supercell (used by the
# `SpinWaveTheory` constructor). Mirrors the plain flatten
# (`reshape_supercell_aux` onto `resize_and_flatten_crystal`), applied to both
# the contracted and physical systems, then rebuilds the mapping so the
# flattened system still carries per-unit `dipole_operators` (needed for the
# Zeeman term in `swt_data`).
flatten_supercell_entangled(sys::System) =
    reshape_entangled(sys, s -> reshape_supercell_aux(s, resize_and_flatten_crystal(s.crystal, s.dims), (1, 1, 1)))

# Apply an ordinary reshaping operation `f` to the *contracted* part of an
# entangled `sys`, treating it as a plain system. The `entanglement` field is
# temporarily detached so `f` (which dispatches on `is_entangled`) takes the
# ordinary path and does not recurse. The reshaped system returned by `f`
# already carries `entanglement = nothing` (built via the positional `System`
# constructor); the caller re-attaches fresh metadata.
function reshape_contracted(sys::System, f)
    ent = sys.entanglement
    sys.entanglement = nothing
    try
        return f(sys)
    finally
        sys.entanglement = ent
    end
end