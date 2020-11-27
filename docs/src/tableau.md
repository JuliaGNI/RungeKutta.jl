# Tableaus

If required, it is straight-forward to create a custom tableau.
The tableau of Heun's method, for example, is defined as follows:
```julia
a = [[0.0  0.0]
     [1.0  0.0]]
b =  [0.5, 0.5]
c =  [0.0, 1.0]
o = 2

tab = Tableau(:heun, o, a, b, c)
```
