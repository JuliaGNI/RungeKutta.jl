# Lobatto III Runge-Kutta Methods

```@setup tableaus
using RungeKutta
using RungeKutta.Tableaus: lobatto_a_coefficients, lobatto_b_coefficients,
                           lobatto_c_coefficients, lobatto_c̄_coefficients,
                           lobatto_d_coefficients, lobatto_e_coefficients,
                           lobatto_f_coefficients, lobatto_g_coefficients
using QuadratureRules: lobatto_legendre_nodes, lobatto_legendre_weights
using Main.SymbolicTableaus: symtype, tableaus_markdown
```

## Lobatto IIIA

```@eval tableaus
tableaus_markdown(2:4) do s
    Tableau(Symbol("LobattoIIIA($s)"), 2s-2, lobatto_a_coefficients(symtype(), s), lobatto_legendre_weights(symtype(), s), lobatto_legendre_nodes(symtype(), s))
end
```

## Lobatto IIIB

```@eval tableaus
tableaus_markdown(2:4) do s
    Tableau(Symbol("LobattoIIIB($s)"), 2s-2, lobatto_b_coefficients(symtype(), s), lobatto_legendre_weights(symtype(), s), lobatto_legendre_nodes(symtype(), s))
end
```

## Lobatto IIIC

```@eval tableaus
tableaus_markdown(2:4) do s
    Tableau(Symbol("LobattoIIIC($s)"), 2s-2, lobatto_c_coefficients(symtype(), s), lobatto_legendre_weights(symtype(), s), lobatto_legendre_nodes(symtype(), s))
end
```

## Lobatto IIIC̄

```@eval tableaus
tableaus_markdown(2:4) do s
    Tableau(Symbol("LobattoIIIC̄($s)"), 2s-2, lobatto_c̄_coefficients(symtype(), s), lobatto_legendre_weights(symtype(), s), lobatto_legendre_nodes(symtype(), s))
end
```

## Lobatto IIID

```@eval tableaus
tableaus_markdown(2:4) do s
    Tableau(Symbol("LobattoIIID($s)"), 2s-2, lobatto_d_coefficients(symtype(), s), lobatto_legendre_weights(symtype(), s), lobatto_legendre_nodes(symtype(), s))
end
```

## Lobatto IIIE

```@eval tableaus
tableaus_markdown(2:4) do s
    Tableau(Symbol("LobattoIIIE($s)"), 2s-2, lobatto_e_coefficients(symtype(), s), lobatto_legendre_weights(symtype(), s), lobatto_legendre_nodes(symtype(), s))
end
```

## Lobatto IIIF

```@eval tableaus
tableaus_markdown(2:4) do s
    Tableau(Symbol("LobattoIIIF($s)"), 2s-2, lobatto_f_coefficients(symtype(), s), lobatto_legendre_weights(symtype(), s), lobatto_legendre_nodes(symtype(), s))
end
```

## Lobatto IIIG

```@eval tableaus
tableaus_markdown(2:4) do s
    Tableau(Symbol("LobattoIIIG($s)"), 2s-2, lobatto_g_coefficients(symtype(), s), lobatto_legendre_weights(symtype(), s), lobatto_legendre_nodes(symtype(), s))
end
```
