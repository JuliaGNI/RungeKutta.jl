# Radau-II Runge-Kutta Methods

```@setup tableaus
using RungeKutta
using RungeKutta.Tableaus: radau_2_coefficients
using QuadratureRules: radau_legendre_nodes, radau_legendre_weights
using Main.SymbolicTableaus: symtype, tableaus_markdown
```

## Radau IIA

```@eval tableaus
tableaus_markdown(2:3) do s
    Tableau(Symbol("RadauIIA($s)"), 2s-1, radau_2_coefficients(symtype(), s), radau_legendre_weights(symtype(), s, Val(:right)), radau_legendre_nodes(symtype(), s, Val(:right)))
end
```

## Radau IIB

```@eval tableaus
tableaus_markdown(2:3) do s
    a = radau_2_coefficients(symtype(),s)
    b = radau_legendre_weights(symtype(), s, Val(:right))
    c = radau_legendre_nodes(symtype(), s, Val(:right))
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau(Symbol("RadauIIB($s)"), 2s-1, ā, b, c)
end
```
