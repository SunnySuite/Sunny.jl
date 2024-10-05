# Maps each Hall number to an Spglib SpacegroupType
const all_spacegroup_types = Spglib.get_spacegroup_type.(1:530)

# Each spacegroup 1..230 is associated with one "conventional" setting (choice
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

# Map a Hall number to the standard setting as a Hall number
function standard_setting_for_hall_number(hall_number)
    Sunny.standard_setting[Sunny.all_spacegroup_types[hall_number].number]
end

# For each Hall number 1..530, define the operation P that transforms to the
# standard setting of the associated spacegroup number. The corresponding Hall
# numbers for these standard settings are listed in `standard_setting`. Data is
# sourced from PyXTal,
# https://github.com/MaterSim/PyXtal/blob/1ba044cace1815d450e476a1fcb2fe8cb5798923/pyxtal/database/HM_Full.csv

const transform_to_standard_setting_repr = [
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

# Returns the affine transformation (P, p) that maps a position (direct lattice
# units) in the setting of hall_number to a position in the standard setting
# (different direct lattice units): xₛ = P x + p.
function transform_to_standard_setting(hall_number)
    op = parse_op(replace(transform_to_standard_setting_repr[hall_number],
                 "a"=>"x", "b"=>"y", "c"=>"z"))
    # A crystal position in global coordinates can be viewed as a contraction of
    # lattice vectors and coefficients, [a, b, c]ᵀ [x, y, z]. Transformation of
    # lattice vectors, [aₛ, bₛ, cₛ] = R [a, b, c] can alternatively be viewed as
    # transformation of positions, [xₛ, yₛ, zₛ] = Rᵀ [x, y, z].
    P = op.R'
    p = op.T
    return SymOp(P, p)
end
