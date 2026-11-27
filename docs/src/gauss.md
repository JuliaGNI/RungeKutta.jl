# Gauss-Legendre Runge-Kutta Methods

```@setup tableaus
using RungeKutta
using RungeKutta.Tableaus: gauss_coefficients
using QuadratureRules: gauss_legendre_nodes, gauss_legendre_weights
using Main.SymbolicTableaus: symtype, tableaus_markdown
```

```@eval tableaus
tableaus_markdown(1:3) do s
    Tableau(Symbol("Gauss($s)"), 2s, gauss_coefficients(symtype(), s), gauss_legendre_weights(symtype(), s), gauss_legendre_nodes(symtype(), s))
end
```
