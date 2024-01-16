using RungeKutta
using Test

import SymPy

"""
symtype()

Return `Sym{T}` for `T` being the underlying type of `Sym(1)`.
"""
symtype() = typeof(SymPy.Sym(1))


include("test_utils.jl")
include("test_tableau.jl")
include("test_tableau_partitioned.jl")
include("test_order_conditions.jl")
include("test_symmetry.jl")
include("test_symplecticity.jl")
include("test_gauss.jl")
include("test_lobatto.jl")
include("test_radau.jl")
include("test_tableaus.jl")
include("test_tableaus_prk.jl")
include("test_list.jl")
