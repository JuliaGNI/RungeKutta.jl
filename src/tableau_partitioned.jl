@doc raw"""
Tableau of a Partitioned Runge-Kutta method
```math
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, v(t_{n} + c_j \Delta t, Q_{n,j}, P_{n,j}) , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i}  \, v(t_{n} + c_j \Delta t, Q_{n,i}, P_{n,i}) , \\
P_{n,i} &= p_{n} + h  \sum \limits_{i=1}^{s} \bar{a}_{ij} \, f(t_{n} + c_j \Delta t, Q_{n,j}, P_{n,j}) , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i}   \, f(t_{n} + c_j \Delta t, Q_{n,i}, P_{n,i}) .
\end{aligned}
```

Parameters:
 * `T`: datatype of coefficient arrays

Fields:
 * `name`: symbolic name of the tableau
 * `o`: order of the method
 * `s`: number of stages
 * `q`: Tableau for `q`
 * `p`: Tableau for `p`
 * `R∞`: stability function at infinity

The actual tableaus are stored in `q` and `p`:
 * `a`: coefficients $a_{ij}$ with $ 1 \le i,j \le s$
 * `b`: weights $b_{i}$  with $ 1 \le i \le s$
 * `c`: nodes $c_{i}$  with $ 1 \le i \le s$

Constructors:
```julia
PartitionedTableau{T}(name, o, q, p)
PartitionedTableau{T}(name, q, p)
PartitionedTableau(name::Symbol, q::Tableau, p::Tableau)
PartitionedTableau(name::Symbol, q::Tableau)
```
"""
struct PartitionedTableau{T, S, RT <: Union{Real,Missing}, RTq, RTp, L} <: AbstractPartitionedTableau{T}
    @TableauHeader

    q::Tableau{T, S, RTq, L}
    p::Tableau{T, S, RTp, L}

    R∞::RT

    function PartitionedTableau{T}(name, o, q, p; R∞=missing) where {T}
        @assert q.s == p.s
        if ismissing(R∞) && !ismissing(q.R∞) && !ismissing(p.R∞)
            if q.R∞ == p.R∞
                R∞ = q.R∞
            end
        end
        new{T, q.s, typeof(R∞), typeof(q.R∞), typeof(p.R∞), q.s * q.s}(name, o, q.s, q, p, R∞)
    end

    function PartitionedTableau{T}(name, q, p; kwargs...) where {T}
        PartitionedTableau{T}(name, min(q.o, p.o), q, p; kwargs...)
    end
end

PartitionedTableau(name::Symbol, q::Tableau{T}, p::Tableau{T}; kwargs...) where {T} = PartitionedTableau{T}(name, q, p; kwargs...)
PartitionedTableau(name::Symbol, q::Tableau) = PartitionedTableau(name, q, q; R∞=q.R∞)
PartitionedTableau(q::Tableau, p::Tableau; kwargs...) = PartitionedTableau(Symbol(q.name, p.name), q, p; kwargs...)
PartitionedTableau(q::Tableau) = PartitionedTableau(q.name, q)

Base.hash(tab::PartitionedTableau, h::UInt) = hash(tab.o, hash(tab.s, hash(tab.q, hash(tab.p, hash(:PartitionedTableau, h)))))

Base.:(==)(tab1::PartitionedTableau, tab2::PartitionedTableau) = (
                                            tab1.o == tab2.o
                                         && tab1.s == tab2.s
                                         && tab1.q == tab2.q
                                         && tab1.p == tab2.p
                                         && ((ismissing(tab1.R∞) && ismissing(tab2.R∞)) || (tab1.R∞ == tab2.R∞)))

Base.isapprox(tab1::PartitionedTableau, tab2::PartitionedTableau; kwargs...) = (
                                            tab1.o == tab2.o
                                         && tab1.s == tab2.s
                                         && ((ismissing(tab1.R∞) && ismissing(tab2.R∞)) || (tab1.R∞ == tab2.R∞))
                                         && isapprox(tab1.q, tab2.q; kwargs...)
                                         && isapprox(tab1.p, tab2.p; kwargs...))
 
Base.isequal(tab1::PartitionedTableau{T1}, tab2::PartitionedTableau{T2}) where {T1,T2} = (tab1 == tab2 && T1 == T2 && tab1.name == tab2.name)

Base.eltype(::PartitionedTableau{T}) where {T} = T

GeometricBase.name(tab::PartitionedTableau) = tab.name
GeometricBase.order(tab::PartitionedTableau) = tab.o
nstages(tab::PartitionedTableau) = tab.s
eachstage(tab::PartitionedTableau) = 1:nstages(tab)

GeometricBase.description(tab::PartitionedTableau) = "$(tab.name) with $(tab.s) $(tab.s == 1 ? "stage" : "stages") and order $(tab.o)"
GeometricBase.reference(tab::PartitionedTableau) = reference(Val(tab.name))

isexplicit(tab::PartitionedTableau) = istril(tab.q.a) && istril(tab.p.a) && all([tab.q.a[i,i] == 0 || tab.p.a[i,i] == 0 for i in 1:tab.s]) && (tab.q.c[1] == 0 || tab.p.c[1] == 0)
isimplicit(tab::PartitionedTableau) = !isexplicit(tab)
isdiagonallyimplicit(tab::PartitionedTableau) = isimplicit(tab) && tab.s != 1 && istril(tab.q.a) && istril(tab.p.a)
isfullyimplicit(tab::PartitionedTableau) = !isexplicit(tab) && !(isdiagonallyimplicit(tab))


"""
```julia
Base.show(io::IO, tab::PartitionedTableau)
```

Pretty-print Partitioned Runge-Kutta tableau.
"""
function Base.show(io::IO, tab::PartitionedTableau)
    print(io, "\nPartitioned Runge-Kutta Tableau $(description(tab)):\n")
    show_coefficients(io, tab.q)
    show_coefficients(io, tab.p)
end
