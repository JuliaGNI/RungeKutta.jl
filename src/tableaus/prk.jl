
"""
Partitioned Gauss-Legendre Runge-Kutta tableau with s stages

```julia
PartitionedTableauGauss(::Type{T}, s)
PartitionedTableauGauss(s) = PartitionedTableauGauss(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This [`PartitionedTableau`](@ref) uses [`TableauGauss`](@ref) for both coefficients `a` and `ā`.
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

This [`PartitionedTableau`](@ref) uses [`TableauLobattoIIIA`](@ref) for `a` and [`TableauLobattoIIIB`](@ref) for `ā`.
"""
function TableauLobattoIIIAIIIB(::Type{T}, s::Int) where {T}
    PartitionedTableau(Symbol("LobattoIIIAIIIB", s), TableauLobattoIIIA(T, s),
        TableauLobattoIIIB(T, s); R∞ = (-1)^(s+1))
end

"""
Tableau for Gauss-Lobatto IIIB-IIIA method with s stages

```julia
TableauLobattoIIIBIIIA(::Type{T}, s)
TableauLobattoIIIBIIIA(s) = TableauLobattoIIIBIIIA(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This [`PartitionedTableau`](@ref) uses [`TableauLobattoIIIB`](@ref) for `a` and [`TableauLobattoIIIA`](@ref) for `ā`.
"""
function TableauLobattoIIIBIIIA(::Type{T}, s::Int) where {T}
    PartitionedTableau(Symbol("LobattoIIIBIIIA", s), TableauLobattoIIIB(T, s),
        TableauLobattoIIIA(T, s); R∞ = (-1)^(s+1))
end

"""
Tableau for Gauss-Lobatto IIIA-IIIĀ method with s stages

```julia
TableauLobattoIIIAIIIĀ(::Type{T}, s)
TableauLobattoIIIAIIIĀ(s) = TableauLobattoIIIAIIIĀ(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This [`PartitionedTableau`](@ref) uses [`TableauLobattoIIIA`](@ref) for `a` and [`TableauLobattoIIIĀ`](@ref) for `ā`.
"""
function TableauLobattoIIIAIIIĀ(::Type{T}, s::Int) where {T}
    PartitionedTableau(Symbol("LobattoIIIAIIIĀ", s), TableauLobattoIIIA(T, s),
        TableauLobattoIIIĀ(T, s); R∞ = (-1)^(s+1))
end

"""
Tableau for Gauss-Lobatto IIIB-IIIB̄ method with s stages

```julia
TableauLobattoIIIBIIIB̄(::Type{T}, s)
TableauLobattoIIIBIIIB̄(s) = TableauLobattoIIIBIIIB̄(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This [`PartitionedTableau`](@ref) uses [`TableauLobattoIIIB`](@ref) for `a` and [`TableauLobattoIIIB̄`](@ref) for `ā`.
"""
function TableauLobattoIIIBIIIB̄(::Type{T}, s::Int) where {T}
    PartitionedTableau(Symbol("LobattoIIIBIIIB̄", s), TableauLobattoIIIB(T, s),
        TableauLobattoIIIB̄(T, s); R∞ = (-1)^(s+1))
end

"""
Tableau for Gauss-Lobatto IIIC-IIIC̄ method with s stages

```julia
TableauLobattoIIICIIIC̄(::Type{T}, s)
TableauLobattoIIICIIIC̄(s) = TableauLobattoIIICIIIC̄(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This [`PartitionedTableau`](@ref) uses [`TableauLobattoIIIC`](@ref) for `a` and [`TableauLobattoIIIC̄`](@ref) for `ā`.
"""
function TableauLobattoIIICIIIC̄(::Type{T}, s::Int) where {T}
    PartitionedTableau(Symbol("LobattoIIICIIIC̄", s), TableauLobattoIIIC(T, s),
        TableauLobattoIIIC̄(T, s); R∞ = (-1)^(s+1))
end

"""
Tableau for Gauss-Lobatto IIIC̄-IIIC method with s stages

```julia
TableauLobattoIIIC̄IIIC(::Type{T}, s)
TableauLobattoIIIC̄IIIC(s) = TableauLobattoIIIC̄IIIC(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This [`PartitionedTableau`](@ref) uses [`TableauLobattoIIIC̄`](@ref) for `a` and [`TableauLobattoIIIC`](@ref) for `ā`.
"""
function TableauLobattoIIIC̄IIIC(::Type{T}, s::Int) where {T}
    PartitionedTableau(Symbol("LobattoIIIC̄IIIC", s), TableauLobattoIIIC̄(T, s),
        TableauLobattoIIIC(T, s); R∞ = (-1)^(s+1))
end

"""
Tableau for Gauss-Lobatto IIID-IIID̄ method with s stages

```julia
TableauLobattoIIIDIIID̄(::Type{T}, s)
TableauLobattoIIIDIIID̄(s) = TableauLobattoIIIDIIID̄(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This [`PartitionedTableau`](@ref) uses [`TableauLobattoIIID`](@ref) for `a` and [`TableauLobattoIIID̄`](@ref) for `ā`.
"""
function TableauLobattoIIIDIIID̄(::Type{T}, s::Int) where {T}
    PartitionedTableau(Symbol("LobattoIIIDIIID̄", s), TableauLobattoIIID(T, s),
        TableauLobattoIIID̄(T, s); R∞ = (-1)^s)
end

"""
Tableau for Gauss-Lobatto IIIE-IIIĒ method with s stages

```julia
TableauLobattoIIIEIIIĒ(::Type{T}, s)
TableauLobattoIIIEIIIĒ(s) = TableauLobattoIIIEIIIĒ(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This [`PartitionedTableau`](@ref) uses [`TableauLobattoIIIE`](@ref) for `a` and [`TableauLobattoIIIĒ`](@ref) for `ā`.
"""
function TableauLobattoIIIEIIIĒ(::Type{T}, s::Int) where {T}
    PartitionedTableau(Symbol("LobattoIIIEIIIĒ", s), TableauLobattoIIIE(T, s),
        TableauLobattoIIIĒ(T, s); R∞ = (-1)^s)
end

"""
Tableau for Gauss-Lobatto IIIF-IIIF̄ method with s stages

```julia
TableauLobattoIIIFIIIF̄(::Type{T}, s)
TableauLobattoIIIFIIIF̄(s) = TableauLobattoIIIFIIIF̄(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This [`PartitionedTableau`](@ref) uses [`TableauLobattoIIIF`](@ref) for `a` and [`TableauLobattoIIIF̄`](@ref) for `ā`.
"""
function TableauLobattoIIIFIIIF̄(::Type{T}, s::Int) where {T}
    PartitionedTableau(Symbol("LobattoIIIFIIIF̄", s), TableauLobattoIIIF(T, s),
        TableauLobattoIIIF̄(T, s); R∞ = (-1)^s)
end

"""
Tableau for Gauss-Lobatto IIIF̄-IIIF method with s stages

```julia
TableauLobattoIIIF̄IIIF(::Type{T}, s)
TableauLobattoIIIF̄IIIF(s) = TableauLobattoIIIF̄IIIF(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This [`PartitionedTableau`](@ref) uses [`TableauLobattoIIIF̄`](@ref) for `a` and [`TableauLobattoIIIF`](@ref) for `ā`.
"""
function TableauLobattoIIIF̄IIIF(::Type{T}, s::Int) where {T}
    PartitionedTableau(Symbol("LobattoIIIF̄IIIF", s), TableauLobattoIIIF̄(T, s),
        TableauLobattoIIIF(T, s); R∞ = (-1)^s)
end

"""
Tableau for Gauss-Lobatto IIIG-IIIḠ method with s stages

```julia
TableauLobattoIIIGIIIḠ(::Type{T}, s)
TableauLobattoIIIGIIIḠ(s) = TableauLobattoIIIGIIIḠ(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

This [`PartitionedTableau`](@ref) uses [`TableauLobattoIIIG`](@ref) for `a` and [`TableauLobattoIIIḠ`](@ref) for `ā`.
"""
function TableauLobattoIIIGIIIḠ(::Type{T}, s::Int) where {T}
    PartitionedTableau(Symbol("LobattoIIIGIIIḠ", s), TableauLobattoIIIG(T, s),
        TableauLobattoIIIḠ(T, s); R∞ = (-1)^s)
end

TableauLobattoIIIAIIIB(s) = TableauLobattoIIIAIIIB(Float64, s)
TableauLobattoIIIBIIIA(s) = TableauLobattoIIIBIIIA(Float64, s)
TableauLobattoIIIAIIIĀ(s) = TableauLobattoIIIAIIIĀ(Float64, s)
TableauLobattoIIIBIIIB̄(s) = TableauLobattoIIIBIIIB̄(Float64, s)
TableauLobattoIIICIIIC̄(s) = TableauLobattoIIICIIIC̄(Float64, s)
TableauLobattoIIIC̄IIIC(s) = TableauLobattoIIIC̄IIIC(Float64, s)
TableauLobattoIIIDIIID̄(s) = TableauLobattoIIIDIIID̄(Float64, s)
TableauLobattoIIIEIIIĒ(s) = TableauLobattoIIIEIIIĒ(Float64, s)
TableauLobattoIIIFIIIF̄(s) = TableauLobattoIIIFIIIF̄(Float64, s)
TableauLobattoIIIF̄IIIF(s) = TableauLobattoIIIF̄IIIF(Float64, s)
TableauLobattoIIIGIIIḠ(s) = TableauLobattoIIIGIIIḠ(Float64, s)
