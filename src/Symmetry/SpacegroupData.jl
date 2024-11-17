# Maps each Hall number to an Spglib SpacegroupType
const all_spacegroup_types = Spglib.get_spacegroup_type.(1:530)

function all_spacegroup_types_for_symbol(sgnum::Int)
    1 <= sgnum <= 230 || error("Spacegroup $sgnum outside range 1..230")
    return filter(all_spacegroup_types) do sgt
        sgt.number == sgnum
    end
end

function all_spacegroup_types_for_symbol(symbol::String)
    symbol = replace(symbol, " "=>"")
    ret = filter(all_spacegroup_types) do sgt
        symbols = (sgt.international_short, sgt.international_full, sgt.hall_symbol)
        symbol in replace.(symbols, " " => "")
    end
    isempty(ret) && error("Unknown spacegroup $symbol")
    return ret
end

function suggestion_to_disambiguate_symbol(symbol)
    sgts = all_spacegroup_types_for_symbol(symbol)
    short_symbols = [sgt.international_short for sgt in sgts]
    full_symbols = [sgt.international_full for sgt in sgts]
    choices = [sgt.choice for sgt in sgts]

    if allunique(short_symbols)
        # Short symbols preferred when unambiguous. For example,
        # spacegroup 230 is better written "Pccm" than "P 2/c 2/c 2/m".
        "Disambiguate with one of: " * repr(short_symbols)
    elseif allunique(full_symbols)
        # Sometimes full symbol is needed. Spacegroup 5 is abbreviated
        # "C2", but requires "C 1 2 1", "A 1 2 1", ... to disambiguate.
        "Disambiguate with one of: " * repr(full_symbols)
    else
        # Origin choice "1" or "2" is sufficient to disambiguate
        @assert choices == ["1", "2"]
        @assert all(sgt -> in(sgt.number, standard_setting_differs_in_spglib), sgts)
        "Disambiguate with additional argument: choice=\"1\" or choice=\"2\""
    end
end

# Get the single SpacegroupType associated with symbol that is valid for
# latvecs. Throw an error if the setting is ambiguous.
function unique_spacegroup_type(symbol, latvecs; choice=nothing)
    if symbol isa Int && isnothing(choice)
        # If only spacegroup number provided, then use ITA standard setting
        sgts = [all_spacegroup_types[standard_setting[symbol]]]
    else
        # Otherwise, look up all possible spacegroups and filter by `choice`
        sgts = all_spacegroup_types_for_symbol(symbol)
        if !isnothing(choice)
            sgts = filter(sgts) do sgt
                sgt.choice == choice
            end
            isempty(sgts) && error("Unknown setting choice \"$choice\" for spacegroup $symbol")
            @assert length(sgts) == 1
        end
    end

    # Validate lattice vectors to give a good error message
    cell = cell_type(latvecs)
    hall_cells = [cell_type(Int(sgt.hall_number)) for sgt in sgts]
    compatible_cells = union(all_compatible_cells.(hall_cells)...)
    if !(cell in compatible_cells)
        expected = join(repr.(unique(hall_cells)), " or ")
        error("Expected $expected cell but found $cell.")
    end

    # Check consistency with shape of latvecs
    sgts = filter(sgts) do sgt
        is_spacegroup_type_consistent(sgt, latvecs)
    end
    if isempty(sgts)
        error("Incompatible $cell cell shape")
    end

    if length(sgts) == 1
        return only(sgts)
    elseif length(sgts) > 1
        @assert isnothing(choice)
        error(suggestion_to_disambiguate_symbol(symbol))
    end
end


# Each spacegroup 1..230 is associated with one "ITA standard" setting (choice
# of axes and origin). Following Volume A of International Tables for
# Crystallography, this setting satisfies:
#
# * Unique axis b and cell choice 1 for monoclinic groups
# * Hexagonal axes for rhombohedral groups
# * Origin choice 2 for certain centrosymmetric spacegroups (orthorhombic,
#   tetragonal or cubic) with two options. With this choice, inversion centers
#   placed at (0,0,0).
#
# The table below provides the conventional Hall number 1..530 for each
# spacegroup number between 1..230. The data is sourced from PyXTal,
# https://github.com/MaterSim/PyXtal/blob/1ba044cace1815d450e476a1fcb2fe8cb5798923/doc/Settings.rst#space-group
const standard_setting = [
    1,   2,   3,   6,   9,   18,  21,  30,  39,  57,  60,  63,  72,  81,  90,
    108, 109, 112, 115, 116, 119, 122, 123, 124, 125, 128, 134, 137, 143, 149,
    155, 161, 164, 170, 173, 176, 182, 185, 191, 197, 203, 209, 212, 215, 218,
    221, 227, 229, 230, 234, 239, 245, 251, 257, 263, 266, 269, 275, 279, 284,
    290, 292, 298, 304, 310, 313, 316, 323, 334, 336, 337, 338, 341, 343, 349,
    350, 351, 352, 353, 354, 355, 356, 357, 358, 360, 362, 363, 365, 366, 367,
    368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382,
    383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397,
    398, 399, 400, 401, 403, 405, 406, 407, 409, 411, 412, 413, 415, 417, 418,
    419, 421, 423, 424, 425, 427, 429, 430, 431, 432, 433, 435, 436, 438, 439,
    440, 441, 442, 443, 444, 446, 447, 448, 449, 450, 452, 454, 455, 456, 457,
    458, 460, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474,
    475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489,
    490, 491, 492, 493, 494, 496, 497, 499, 500, 501, 502, 503, 504, 505, 506,
    507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 519, 520, 522, 523,
    524, 526, 528, 529, 530
]

# The Spglib convention for standard setting differs slightly from the ITA one
# above. For each spacegroup number, Spglib selects the first available setting
# in the order of Hall numbers. This corresponds to origin choice 1 rather than
# origin choice 2 for certain centrosymmetric groups with two choices. For
# example, Hall number 525 (instead of 526) will be chosen for the space group
# 227.
const standard_setting_for_spglib = map(1:230) do sgnum
    first(all_spacegroup_types_for_symbol(sgnum)).hall_number
end

# Spacegroup numbers [48, 50, 59, ..., 227, 228] where the Spglib standard
# setting differs from the ITA one.
const standard_setting_differs_in_spglib = filter(1:230) do sgnum
    standard_setting[sgnum] != standard_setting_for_spglib[sgnum]
end


# Map a Hall number to the standard setting as a Hall number
function standard_setting_for_hall_number(hall_number)
    Sunny.standard_setting[Sunny.all_spacegroup_types[hall_number].number]
end

# For each Hall number 1..530, a string representation of the affine
# transformation that maps lattice vectors to the ITA standard setting of the
# associated spacegroup number. The Hall numbers of these standard settings are
# listed in `standard_setting`. Data is sourced from PyXTal,
# https://github.com/MaterSim/PyXtal/blob/1ba044cace1815d450e476a1fcb2fe8cb5798923/pyxtal/database/HM_Full.csv
const mapping_to_standard_setting_repr = [
    "a,b,c", "a,b,c", "a,b,c", "c,a,b", "b,c,a", "a,b,c", "c,a,b", "b,c,a",
    "a,b,c", "c,-b,a", "c,b,-a-c", "c,a,b", "a,c,-b", "-a-c,c,b", "b,c,a",
    "-b,a,c", "b,-a-c,c", "a,b,c", "c,a,b", "b,c,a", "a,b,c", "-a-c,b,a",
    "c,-b,a", "c,a,b", "a,-a-c,b", "a,c,-b", "b,c,a", "b,a,-a-c", "-b,a,c",
    "a,b,c", "c,-b,a", "c,b,-a-c", "c,a,b", "a,c,-b", "-a-c,c,b", "b,c,a",
    "-b,a,c", "b,-a-c,c", "a,b,c", "-a-c,b,a", "c,b,-a-c", "c,-b,a",
    "a,-b,-a-c", "-a-c,-b,c", "c,a,b", "a,-a-c,b", "-a-c,c,b", "a,c,-b",
    "-a-c,a,-b", "c,-a-c,-b", "b,c,a", "b,a,-a-c", "b,-a-c,c", "-b,a,c",
    "-b,-a-c,a", "-b,c,-a-c", "a,b,c", "c,a,b", "b,c,a", "a,b,c", "c,a,b",
    "b,c,a", "a,b,c", "c,-b,a", "c,b,-a-c", "c,a,b", "a,c,-b", "-a-c,c,b",
    "b,c,a", "-b,a,c", "b,-a-c,c", "a,b,c", "-a-c,b,a", "c,-b,a", "c,a,b",
    "a,-a-c,b", "a,c,-b", "b,c,a", "b,a,-a-c", "-b,a,c", "a,b,c", "-a-c,b,a",
    "c,-b,a", "c,a,b", "a,-a-c,b", "-a-c,c,b", "b,c,a", "b,a,-a-c", "b,-a-c,c",
    "a,b,c", "-a-c,b,a", "c,b,-a-c", "c,-b,a", "a,-b,-a-c", "-a-c,-b,c",
    "c,a,b", "a,-a-c,b", "-a-c,c,b", "a,c,-b", "-a-c,a,-b", "c,-a-c,-b",
    "b,c,a", "b,a,-a-c", "b,-a-c,c", "-b,a,c", "-b,-a-c,a", "-b,c,-a-c",
    "a,b,c", "a,b,c", "c,a,b", "b,c,a", "a,b,c", "c,a,b", "b,c,a", "a,b,c",
    "a,b,c", "c,a,b", "b,c,a", "a,b,c", "c,a,b", "b,c,a", "a,b,c", "a,b,c",
    "a,b,c", "a,b,c", "c,a,b", "b,c,a", "a,b,c", "b,a,-c", "c,a,b", "-c,b,a",
    "b,c,a", "a,-c,b", "a,b,c", "c,a,b", "b,c,a", "a,b,c", "b,a,-c", "c,a,b",
    "-c,b,a", "b,c,a", "a,-c,b", "a,b,c", "b,a,-c", "c,a,b", "-c,b,a", "b,c,a",
    "a,-c,b", "a,b,c", "b,a,-c", "c,a,b", "-c,b,a", "b,c,a", "a,-c,b", "a,b,c",
    "b,a,-c", "c,a,b", "-c,b,a", "b,c,a", "a,-c,b", "a,b,c", "c,a,b", "b,c,a",
    "a,b,c", "b,a,-c", "c,a,b", "-c,b,a", "b,c,a", "a,-c,b", "a,b,c", "c,a,b",
    "b,c,a", "a,b,c", "c,a,b", "b,c,a", "a,b,c", "b,a,-c", "c,a,b", "-c,b,a",
    "b,c,a", "a,-c,b", "a,b,c", "c,a,b", "b,c,a", "a,b,c", "b,a,-c", "c,a,b",
    "-c,b,a", "b,c,a", "a,-c,b", "a,b,c", "b,a,-c", "c,a,b", "-c,b,a", "b,c,a",
    "a,-c,b", "a,b,c", "b,a,-c", "c,a,b", "-c,b,a", "b,c,a", "a,-c,b", "a,b,c",
    "b,a,-c", "c,a,b", "-c,b,a", "b,c,a", "a,-c,b", "a,b,c", "c,a,b", "b,c,a",
    "a,b,c", "c,a,b", "b,c,a", "a,b,c", "c,a,b", "b,c,a", "a,b,c", "c,a,b",
    "b,c,a", "a,b,c", "b,a,-c", "c,a,b", "-c,b,a", "b,c,a", "a,-c,b", "a,b,c",
    "a+1/4,b+1/4,c+1/4", "a,b,c", "a,b,c", "c,a,b", "b,c,a", "a+1/4,b+1/4,c",
    "a,b,c", "c+1/4,a+1/4,b", "c,a,b", "b+1/4,c+1/4,a", "b,c,a", "a,b,c",
    "b,a,-c", "c,a,b", "-c,b,a", "b,c,a", "a,-c,b", "a,b,c", "b,a,-c", "c,a,b",
    "-c,b,a", "b,c,a", "a,-c,b", "a,b,c", "b,a,-c", "c,a,b", "-c,b,a", "b,c,a",
    "a,-c,b", "a,b,c", "b,a,-c", "c,a,b", "-c,b,a", "b,c,a", "a,-c,b", "a,b,c",
    "c,a,b", "b,c,a", "a,b,c", "c,a,b", "b,c,a", "a,b,c", "b,a,-c", "c,a,b",
    "-c,b,a", "b,c,a", "a,-c,b", "a,b,c", "c,a,b", "b,c,a", "a+1/4,b+1/4,c",
    "a,b,c", "c+1/4,a+1/4,b", "c,a,b", "b+1/4,c+1/4,a", "b,c,a", "a,b,c",
    "b,a,-c", "c,a,b", "-c,b,a", "b,c,a", "a,-c,b", "a,b,c", "b,a,-c", "a,b,c",
    "b,a,-c", "c,a,b", "-c,b,a", "b,c,a", "a,-c,b", "a,b,c", "b,a,-c", "c,a,b",
    "-c,b,a", "b,c,a", "a,-c,b", "a,b,c", "b,a,-c", "c,a,b", "-c,b,a", "b,c,a",
    "a,-c,b", "a,b,c", "c,a,b", "b,c,a", "a,b,c", "c,a,b", "b,c,a", "a,b,c",
    "b,a,-c", "c,a,b", "-c,b,a", "b,c,a", "a,-c,b", "a,b+1/4,c+1/4", "a,b,c",
    "b,a+1/4,-c+1/4", "b,a,-c", "c,a+1/4,b+1/4", "c,a,b", "-c,b+1/4,a+1/4",
    "-c,b,a", "b,c+1/4,a+1/4", "b,c,a", "a,-c+1/4,b+1/4", "a,-c,b", "a,b,c",
    "a-1/8,b-1/8,c-1/8", "a,b,c", "a,b,c", "a,b,c", "c,a,b", "b,c,a", "a,b,c",
    "b,a,-c", "a,b,c", "b,a,-c", "c,a,b", "-c,b,a", "b,c,a", "a,-c,b", "a,b,c",
    "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c",
    "a,b,c", "a-1/4,b+1/4,c", "a,b,c", "a-1/4,b-1/4,c-1/4", "a,b,c", "a,b,c",
    "a,b-1/4,c-1/8", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c",
    "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c",
    "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c",
    "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c",
    "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c",
    "a-1/4,b-1/4,c", "a,b,c", "a-1/4,b-1/4,c-1/4", "a,b,c", "a,b,c", "a,b,c",
    "a-1/4,b+1/4,c", "a,b,c", "a-1/4,b+1/4,c", "a,b,c", "a,b,c", "a,b,c",
    "a-1/4,b+1/4,c-1/4", "a,b,c", "a-1/4,b+1/4,c-1/4", "a,b,c", "a,b,c",
    "a,b,c", "a-1/4,b+1/4,c-1/4", "a,b,c", "a-1/4,b+1/4,c-1/4", "a,b,c",
    "a,b,c", "a,b,c", "a,b+1/4,c-1/8", "a,b,c", "a,b+1/4,c-1/8", "a,b,c",
    "a,b,c", "a,b,c", "a,b,c", "a,b,c",
    "2/3a+1/3b+1/3c,-1/3a+1/3b+1/3c,-1/3a-2/3b+1/3c", "a,b,c", "a,b,c",
    "2/3a+1/3b+1/3c,-1/3a+1/3b+1/3c,-1/3a-2/3b+1/3c", "a,b,c", "a,b,c", "a,b,c",
    "a,b,c", "a,b,c", "a,b,c", "a,b,c",
    "2/3a+1/3b+1/3c,-1/3a+1/3b+1/3c,-1/3a-2/3b+1/3c", "a,b,c", "a,b,c", "a,b,c",
    "a,b,c", "a,b,c", "2/3a+1/3b+1/3c,-1/3a+1/3b+1/3c,-1/3a-2/3b+1/3c", "a,b,c",
    "2/3a+1/3b+1/3c,-1/3a+1/3b+1/3c,-1/3a-2/3b+1/3c", "a,b,c", "a,b,c", "a,b,c",
    "a,b,c", "a,b,c", "2/3a+1/3b+1/3c,-1/3a+1/3b+1/3c,-1/3a-2/3b+1/3c", "a,b,c",
    "2/3a+1/3b+1/3c,-1/3a+1/3b+1/3c,-1/3a-2/3b+1/3c", "a,b,c", "a,b,c", "a,b,c",
    "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c",
    "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c",
    "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c",
    "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a-1/4,b-1/4,c-1/4",
    "a,b,c", "a,b,c", "a-1/8,b-1/8,c-1/8", "a,b,c", "a,b,c", "a,b,c", "a,b,c",
    "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c",
    "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c", "a,b,c",
    "a-1/4,b-1/4,c-1/4", "a,b,c", "a,b,c", "a-1/4,b-1/4,c-1/4", "a,b,c",
    "a,b,c", "a,b,c", "a-1/8,b-1/8,c-1/8", "a,b,c", "a-3/8,b-3/8,c-3/8",
    "a,b,c", "a,b,c", "a,b,c"
]

# Returns the affine transformation `setting` that maps a position x in the
# setting of hall_number to a position xₛ in the ITA standard setting for the
# spacegroup: xₛ = transform(setting, x).
function mapping_to_standard_setting(hall_number)
    # These steps were determined empirically by comparison with the Bilbao
    # server, which lists setting strings and P matrices simultaneously.
    op = parse_op(replace(mapping_to_standard_setting_repr[hall_number],
                          "a"=>"x", "b"=>"y", "c"=>"z"))
    P = op.R'
    p = op.T
    return SymOp(P, p)
end

# Given a spacegroup number and a table of symops, try to infer the affine map
# that transforms to the ITA standard setting.
function hall_number_from_symops(sgnum, symops)
    sgts = filter(all_spacegroup_types_for_symbol(sgnum)) do sgt
        Rs, Ts = Spglib.get_symmetry_from_database(sgt.hall_number)
        SymOp.(Rs, Ts) ≈ symops
    end

    if isempty(sgts)
        # Cannot be matched to any of the Hall number settings
        return nothing
    else
        return Int(only(sgts).hall_number)
    end
end

# Returns the affine map from an arbitrary Spglib-inferred dataset to the ITA
# standard setting. Beware that this function may lead to some strange
# behaviors. For example, suppose one initializes diamond-cubic positions
# (spacegroup 227) under the ITA standard setting (choice=2). When Spglib infers
# the setting (assuming choice=1), this function will lead to an origin shift of
# [0.5, 0, 0]; this, in turn, will effectively swap the Wyckoffs 8a and 8b.
function mapping_to_standard_setting_from_spglib_dataset(d::Spglib.Dataset)
    # Map from an arbitrary setting (inferred by Spglib) to the Spglib standard
    # setting. Commonly the elements will be a simple fractions, so search for
    # this case up to 15 digits.
    spglib_from_any = SymOp(
        rationalize.(d.transformation_matrix; tol=1e-15),
        rationalize.(d.origin_shift; tol=1e-15),
    )

    # Map from the Spglib standard setting to the ITA one
    hall_spglib = standard_setting_for_spglib[d.spacegroup_number]
    std_from_spglib = mapping_to_standard_setting(hall_spglib)

    # Return the composed mapping
    return std_from_spglib * spglib_from_any
end


# The lattice system that is expected for a given Hall number. See:
# http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37
function cell_type(hall_number::Int)
    if 1 <= hall_number <= 2
        triclinic
    elseif 3 <= hall_number <= 107
        monoclinic
    elseif 108 <= hall_number <= 348
        orthorhombic
    elseif 349 <= hall_number <= 429
        tetragonal
    elseif 430 <= hall_number <= 461
        # The trigonal space groups require either rhombohedral or hexagonal
        # cells. The Hall numbers below have "choice" R.
        hall_number in [434, 437, 445, 451, 453, 459, 461] ? rhombohedral : hexagonal
    elseif 462 <= hall_number <= 488
        hexagonal
    elseif 489 <= hall_number <= 530
        cubic
    else
        error("Invalid Hall number $hall_number. Allowed range is 1..530")
    end
end

function is_trigonal_symmetry(hall_number::Int)
    return 430 <= hall_number <= 461
end

# Centering type for Hall number. Possible values: 'P' (simple), 'C', 'A'
# (Base), 'I' (Body), 'F' (Face), 'R' (Rhombohedral).
function centering_symbol(hall_number::Int)
    first(all_spacegroup_types[hall_number].international_short)
end

# Primitive basis in multiples of the lattice vectors for the ITA standard
# setting.
const standard_primitive_basis = Dict(
    'P' => SA[1 0 0; 0 1 0; 0 0 1],
    'C' => SA[1/2 1/2 0; -1/2 1/2 0; 0 0 1],
    'A' => SA[1 0 0; 0 1/2 -1/2; 0 1/2 1/2],
    'I' => SA[-1/2 1/2 1/2; 1/2 -1/2 1/2; 1/2 1/2 -1/2],
    'F' => SA[0 1/2 1/2; 1/2 0 1/2; 1/2 1/2 0],
    'R' => SA[2/3 -1/3 -1/3; 1/3 1/3 -2/3; 1/3 1/3 1/3],
)

struct Spacegroup
    symops    :: Vector{SymOp} # Symmetry operations
    label     :: String        # Description of space group
    number    :: Int           # International spacegroup number (1..230)

    # Mapping to ITA standard setting, xₛ = transform(setting, x). The variable
    # x denotes position in direct lattice units for the current crystal
    # (arbitrary setting), and xₛ denotes position in direct lattice units for
    # the standard setting of spacegroup `number`. This becomes xₛ = P x + p
    # where P ≡ setting.R and p ≡ setting.T. Lattice vectors therefore transform
    # as [aₛ bₛ cₛ] = [a b c] P⁻¹ to ensure that positions in global coordinates
    # [a b c] x = [aₛ bₛ cₛ] xₛ are invariant to a change of setting, up to
    # overall translation.
    setting :: SymOp
end

function Spacegroup(hall_number::Int)
    @assert 1 <= hall_number <= 530
    symops = SymOp.(Spglib.get_symmetry_from_database(hall_number)...)
    label = spacegroup_label(hall_number)
    number = Int(all_spacegroup_types[hall_number].number)
    setting = mapping_to_standard_setting(hall_number)
    return Spacegroup(symops, label, number, setting)
end
