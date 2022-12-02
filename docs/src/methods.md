# Methods with Tabulated Coefficients

RungeKutta.jl provides tabulated coefficients for various explicit and implicit (both diagonally and fully implicit) Runge-Kutta methods, in particular

| Function and Aliases                                           | Stages | Order |
|:---------------------------------------------------------------|:-------|:------|
| **Explicit Methods**                                           |        |       |
| [`TableauExplicitEuler`](@ref), [`TableauForwardEuler`](@ref)  | 1      | 1     |
| [`TableauExplicitMidpoint`](@ref)                              | 2      | 2     |
| [`TableauHeun2`](@ref)                                         | 2      | 2     |
| [`TableauHeun3`](@ref)                                         | 3      | 3     |
| [`TableauRalston2`](@ref)                                      | 2      | 2     |
| [`TableauRalston3`](@ref)                                      | 3      | 3     |
| [`TableauRunge`](@ref), [`TableauRunge2`](@ref)                | 2      | 2     |
| [`TableauKutta`](@ref), [`TableauKutta3`](@ref)                | 3      | 3     |
| [`TableauRK21`](@ref)                                          | 2      | 2     |
| [`TableauRK22`](@ref)                                          | 2      | 2     |
| [`TableauRK31`](@ref)                                          | 3      | 3     |
| [`TableauRK32`](@ref)                                          | 3      | 3     |
| [`TableauRK4`](@ref)                                           | 6      | 5     |
| [`TableauRK41`](@ref)                                          | 4      | 4     |
| [`TableauRK42`](@ref)                                          | 4      | 4     |
| [`TableauRK416`](@ref)                                         | 4      | 4     |
| [`TableauRK438`](@ref)                                         | 4      | 4     |
| [`TableauRK5`](@ref)                                           | 6      | 5     |
| [`TableauSSPRK3`](@ref)                                        | 3      | 3     |
| **Diagonally Implicit Methods**                                |        |       |
| [`TableauCrankNicolson`](@ref)                                 | 2      | 2     |
| [`TableauKraaijevangerSpijker`](@ref)                          | 2      | 2     |
| [`TableauQinZhang`](@ref)                                      | 2      | 2     |
| [`TableauCrouzeix`](@ref)                                      | 2      | 3     |
| **Fully Implicit Methods**                                     |        |       |
| [`TableauImplicitEuler`](@ref), [`TableauBackwardEuler`](@ref) | 1      | 1     |
| [`TableauImplicitMidpoint`](@ref)                              | 2      | 2     |
| [`TableauSRK3`](@ref)                                          | 3      | 4     |

All constructors take an optional type argument, defaulting to `Float64`.


## Gauß, Lobatto and Radau Methods

The coefficients of the Gauß, Lobatto and Radau methods are computed on-the-fly by the following constructors:

| Function                            | Method                      | Order   |
|:------------------------------------|:----------------------------|:--------|
| [`TableauGauss(s)`](@ref)           | Gauß with s stages          | $2s$    |
| [`TableauLobattoIII(s)`](@ref)      | Lobatto III with s stages   | $2s-2$  |
| [`TableauLobattoIIIA(s)`](@ref)     | Lobatto IIIA with s stages  | $2s-2$  |
| [`TableauLobattoIIIB(s)`](@ref)     | Lobatto IIIB with s stages  | $2s-2$  |
| [`TableauLobattoIIIC(s)`](@ref)     | Lobatto IIIC with s stages  | $2s-2$  |
| [`TableauLobattoIIID(s)`](@ref)     | Lobatto IIID with s stages  | $2s-2$  |
| [`TableauLobattoIIIE(s)`](@ref)     | Lobatto IIIE with s stages  | $2s-2$  |
| [`TableauLobattoIIIF(s)`](@ref)     | Lobatto IIIF with s stages  | $2s$    |
| [`TableauLobattoIIIG(s)`](@ref)     | Lobatto IIIG with s stages  | $2s$    |
| [`TableauRadauIA(s)`](@ref)         | Radau IA with s stages      | $2s-1$  |
| [`TableauRadauIB(s)`](@ref)         | Radau IB with s stages      | $2s-1$  |
| [`TableauRadauIIA(s)`](@ref)        | Radau IIA with s stages     | $2s-1$  |
| [`TableauRadauIIB(s)`](@ref)        | Radau IIB with s stages     | $2s-1$  |

The argument `s` refers to the number of stages ($s \ge 1$ for Gauß and $s \ge 2$ for all other methods). The type specifier can also be ommitted.
For each method, a second constructor exists, where the first argument specifies the number type of the coefficients and the second argument the number of stages. Internally, all coefficients are computed using `BigFloat` and then converted to the requested number type, defaulting to `Float64`.
