using RungeKutta

# Prints the full-precision `a`, `b` and `c` coefficients of every exported
# tableau constructor whose only argument is the element type, for Float64 and
# BigFloat. Aliases are left out, and so is `TableauSSPRK2`, which returns the
# coefficients of `TableauHeun2`. Run the script on two commits and diff the
# output: identical output means that the coefficients are bit-identical.

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
