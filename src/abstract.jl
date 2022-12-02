
abstract type AbstractTableau{T} end
abstract type AbstractPartitionedTableau{T} <: AbstractTableau{T} end

@define TableauHeader begin
    name::Symbol
    o::Int
    s::Int
end

name(::AbstractTableau) = missing
order(::AbstractTableau) = missing
nstages(::AbstractTableau) = missing
