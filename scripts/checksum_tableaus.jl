using RungeKutta

# Verifies that removing the stale `big` import in src/Tableaus.jl does not
# change the numerical value of any tableau's coefficients. `@big` controls
# the precision of tableau coefficients at macro-expansion time, so this
# prints the full-precision `a`, `b`, `c` coefficients of every nullary
# tableau constructor (those taking only an optional element type) for both
# Float64 and BigFloat. Run on the base commit and on the branch and diff the
# output: identical output means the coefficients are bit-identical.

const CONSTRUCTORS = [
    TableauCrankNicolson, TableauKraaijevangerSpijker, TableauQinZhang, TableauCrouzeix,
    TableauExplicitEuler, TableauExplicitMidpoint, TableauHeun2, TableauHeun3,
    TableauRalston2, TableauRalston3, TableauRunge, TableauKutta, TableauRK31,
    TableauRK416, TableauRK42, TableauRK438, TableauRK5, TableauSSPRK3,
    TableauImplicitEuler, TableauImplicitMidpoint, TableauIRK3, TableauSRK3
]

for f in CONSTRUCTORS
    for T in (Float64, BigFloat)
        tab = f(T)
        println(nameof(f), "{", T, "}")
        println("  a = ", repr(tab.a))
        println("  b = ", repr(tab.b))
        println("  c = ", repr(tab.c))
    end
end
