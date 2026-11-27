using RungeKutta
using LinearAlgebra

# Checks that `TableauRK42` has order four (Butcher 2016, p. 102), in two independent ways:
#
# 1. the conditions of the rooted trees up to order four hold in exact rational arithmetic, and
#    the row sums of `a` equal `c`;
# 2. the error at t = 1 on y' = y, y(0) = 1, falls with the step size h like h⁴.
#
# The same checks run on the tableau with a[4,3] = 1 as a negative control: it must fail both,
# with an observed order of about one.

function rooted_tree_residuals(a, b, c)
    [b' * ones(eltype(b), length(b)) - 1,
        b' * c - 1 // 2,
        b' * c .^ 2 - 1 // 3,
        b' * (a * c) - 1 // 6,
        b' * c .^ 3 - 1 // 4,
        b' * (c .* (a * c)) - 1 // 8,
        b' * (a * c .^ 2) - 1 // 12,
        b' * (a * (a * c)) - 1 // 24]
end

function global_error(a, b, h)
    s = length(b)
    y = 1.0
    for _ in 1:round(Int, 1 / h)
        k = zeros(s)
        for i in 1:s
            k[i] = y + h * sum(a[i, j] * k[j] for j in 1:(i - 1); init = 0.0)
        end
        y += h * dot(b, k)
    end
    abs(y - exp(1.0))
end

function observed_orders(a, b)
    errs = [global_error(Float64.(a), Float64.(b), h) for h in (1 / 20, 1 / 40, 1 / 80)]
    log2.(errs[1:(end - 1)] ./ errs[2:end])
end

# The coefficients of TableauRK42 are small rationals, which `rationalize` recovers exactly.
tab = TableauRK42(Float64)
a = Matrix(rationalize.(BigInt, tab.a))
b = rationalize.(BigInt, collect(tab.b))
c = rationalize.(BigInt, collect(tab.c))

a_old = copy(a)
a_old[4, 3] = 1

for (label, A) in (("TableauRK42", a), ("negative control, a[4,3] = 1", a_old))
    r = rooted_tree_residuals(A, b, c)
    println(label)
    println("  row sums of a equal c:     ", vec(sum(A; dims = 2)) == c)
    println("  rooted-tree residuals:     ", r)
    println("  observed order (h → h/2):  ", round.(observed_orders(A, b); digits = 2))
end
