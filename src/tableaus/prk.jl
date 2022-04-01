
"""
Partitioned Gauss-Legendre Runge-Kutta tableau with s stages

```julia
PartitionedTableauGauss(::Type{T}, s)
PartitionedTableauGauss(s) = PartitionedTableauGauss(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This partitioned tableau uses [`TableauGauss`](@ref) for both coefficients `a` and `ā`.
"""
function PartitionedTableauGauss(::Type{T}, s::Int) where {T}
    PartitionedTableau(Symbol("PGauss", s), TableauGauss(T, s))
end

PartitionedTableauGauss(s) = PartitionedTableauGauss(Float64, s)


"""
Partitioned Gauss-Lobatto IIIA-IIIB tableau with s stages

```julia
TableauLobattoIIIAIIIB(::Type{T}, s)
TableauLobattoIIIAIIIB(s) = TableauLobattoIIIAIIIB(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This partitioned tableau uses [`TableauLobattoIIIA`](@ref) for `a` and [`TableauLobattoIIIB`](@ref) for `ā`.
"""
function TableauLobattoIIIAIIIB(s)
    PartitionedTableau(Symbol("LobattoIIIAIIIB", s), TableauLobattoIIIA(s), TableauLobattoIIIB(s))
end


"""
Tableau for Gauss-Lobatto IIIB-IIIA method with s stages

```julia
TableauLobattoIIIBIIIA(::Type{T}, s)
TableauLobattoIIIBIIIA(s) = TableauLobattoIIIBIIIA(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This partitioned tableau uses [`TableauLobattoIIIB`](@ref) for `a` and [`TableauLobattoIIIA`](@ref) for `ā`.
"""
function TableauLobattoIIIBIIIA(s)
    PartitionedTableau(Symbol("LobattoIIIBIIIA", s), TableauLobattoIIIB(s), TableauLobattoIIIA(s))
end


"""
Tableau for Gauss-Lobatto IIIC-IIIC̄ method with s stages

```julia
TableauLobattoIIICIIIC̄(::Type{T}, s)
TableauLobattoIIICIIIC̄(s) = TableauLobattoIIICIIIC̄(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This partitioned tableau uses [`TableauLobattoIIIC`](@ref) for `a` and [`TableauLobattoIIIC̄`](@ref) for `ā`.
"""
function TableauLobattoIIICIIIC̄(s)
    PartitionedTableau(Symbol("LobattoIIICIIIC̄", s), TableauLobattoIIIC(s), TableauLobattoIIIC̄(s))
end


"""
Tableau for Gauss-Lobatto IIIC̄-IIIC method with s stages

```julia
TableauLobattoIIIC̄IIIC(::Type{T}, s)
TableauLobattoIIIC̄IIIC(s) = TableauLobattoIIIC̄IIIC(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This partitioned tableau uses [`TableauLobattoIIIC̄`](@ref) for `a` and [`TableauLobattoIIIC`](@ref) for `ā`.
"""
function TableauLobattoIIIC̄IIIC(s)
    PartitionedTableau(Symbol("LobattoIIIC̄IIIC", s), TableauLobattoIIIC̄(s), TableauLobattoIIIC(s))
end
