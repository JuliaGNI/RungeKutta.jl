
abstract type AbstractTableau{T} end
abstract type AbstractPartitionedTableau{T} <: AbstractTableau{T} end

@define TableauHeader begin
    name::Symbol
    o::Int
    s::Int
end

GeometricBase.name(::AbstractTableau) = missing
GeometricBase.order(::AbstractTableau) = missing
nstages(::AbstractTableau) = missing
