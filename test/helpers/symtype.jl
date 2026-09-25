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
function Base.convert(::Type{T}, x::SymPyPythonCall.Sym) where {T <: AbstractFloat}
    T(SymPyPythonCall.N(x))
end
