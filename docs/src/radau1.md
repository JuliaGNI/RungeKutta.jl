# Radau-I Runge-Kutta Methods

```@setup tableaus
using RungeKutta
using RungeKutta.Tableaus: radau_1_coefficients
using QuadratureRules: radau_legendre_nodes, radau_legendre_weights
using Main.SymbolicTableaus: symtype, tableaus_markdown
```

## Radau IA

```@eval tableaus
tableaus_markdown(2:3) do s
    Tableau(Symbol("RadauIA($s)"), 2s-1, radau_1_coefficients(symtype(), s), radau_legendre_weights(symtype(), s, Val(:left)), radau_legendre_nodes(symtype(), s, Val(:left)))
end
```

## Radau IB

```@eval tableaus
tableaus_markdown(2:3) do s
    a = radau_1_coefficients(symtype(),s)
    b = radau_legendre_weights(symtype(), s, Val(:left))
    c = radau_legendre_nodes(symtype(), s, Val(:left))
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau(Symbol("RadauIB($s)"), 2s-1, ā, b, c)
end
```
