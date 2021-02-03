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

The actual tableaus are stored in `q` and `p`:
 * `a`: coefficients $a_{ij}$ with $ 1 \le i,j \le s$
 * `b`: weights $b_{i}$  with $ 1 \le i \le s$
 * `c`: nodes $c_{i}$  with $ 1 \le i \le s$

Constructors:
```julia
PartitionedTableau{T}(name, o, s, q, p)
PartitionedTableau{T}(name, q, p)
PartitionedTableau(name::Symbol, q::Tableau, p::Tableau)
PartitionedTableau(name::Symbol, q::Tableau)
```
"""
struct PartitionedTableau{T}
    name::Symbol
    o::Int
    s::Int

    q::Tableau{T}
    p::Tableau{T}

    function PartitionedTableau{T}(name, o, s, q, p) where {T}
        @assert s == q.s == p.s
        new(name, o, s, q, p)
    end

    function PartitionedTableau{T}(name, q, p) where {T}
        new(name, min(q.o, p.o), q.s, q, p)
    end
end

PartitionedTableau(name::Symbol, q::Tableau{T}, p::Tableau{T}) where {T} = PartitionedTableau{T}(name, q, p)
PartitionedTableau(name::Symbol, q::Tableau) = PartitionedTableau(name, q, q)

Base.hash(tab::PartitionedTableau, h::UInt) = hash(tab.o, hash(tab.s, hash(tab.q, hash(tab.p, hash(:PartitionedTableau, h)))))

Base.:(==)(tab1::PartitionedTableau, tab2::PartitionedTableau) = (
                                            tab1.o == tab2.o
                                         && tab1.s == tab2.s
                                         && tab1.q == tab2.q
                                         && tab1.p == tab2.p)

Base.isapprox(tab1::PartitionedTableau, tab2::PartitionedTableau; kwargs...) = (
                                            tab1.o == tab2.o
                                         && tab1.s == tab2.s
                                         && isapprox(tab1.q, tab2.q; kwargs...)
                                         && isapprox(tab1.p, tab2.p; kwargs...))
 
Base.isequal(tab1::PartitionedTableau{T1}, tab2::PartitionedTableau{T2}) where {T1,T2} = (tab1 == tab2 && T1 == T2 && tab1.name == tab2.name)

Base.eltype(::PartitionedTableau{T}) where {T} = T

name(tab::PartitionedTableau) = tab.name
order(tab::PartitionedTableau) = tab.o
nstages(tab::PartitionedTableau) = tab.s
eachstage(tab::PartitionedTableau) = 1:tab.s

isexplicit(tab::PartitionedTableau) = istril(tab.q.a) && istril(tab.p.a) && all([tab.q.a[i,i] == 0 || tab.p.a[i,i] == 0 for i in 1:tab.s]) && (tab.q.c[1] == 0 || tab.p.c[1] == 0)
isimplicit(tab::PartitionedTableau) = !isexplicit(tab)
isdiagnonallyimplicit(tab::PartitionedTableau) = isimplicit(tab) && tab.s != 1 && istril(tab.q.a) && istril(tab.p.a)
isfullyimplicit(tab::PartitionedTableau) = !isexplicit(tab) && !(isdiagnonallyimplicit(tab))
