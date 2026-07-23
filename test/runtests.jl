using RungeKutta
using Test

import SymPyPythonCall

"""
symtype()

Return `Sym{T}` for `T` being the underlying type of `Sym(1)`.
"""
symtype() = typeof(SymPyPythonCall.Sym(1))

# Convert a symbolic value to a float by evaluating it. SymPyPythonCall's
# `convert(Float64, ::Sym)` uses `pyconvert`, which rejects unevaluated
# expressions (e.g. `1/2 - sqrt(3)/6`); the `N` path evaluates them. This is
# needed for the `symtype() ≈ Float64` comparisons below, whose `isapprox`
# goes through `convert(Float64, ::Sym)` inside `LinearAlgebra.norm`.
Base.convert(::Type{T}, x::SymPyPythonCall.Sym) where {T<:AbstractFloat} = T(SymPyPythonCall.N(x))


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
