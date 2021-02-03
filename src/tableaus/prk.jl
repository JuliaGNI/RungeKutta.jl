
"Gauss-Legendre Runge-Kutta"
function PartitionedTableauGauss(s::Int)
    PartitionedTableau(Symbol("PGauss", s), TableauGauss(s))
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with s stages"
function TableauLobattoIIIAIIIB(s)
    PartitionedTableau(Symbol("LobattoIIIAIIIB", s), TableauLobattoIIIA(s), TableauLobattoIIIB(s))
end

"Tableau for Gauss-Lobatto IIIB-IIIA method with s stages"
function TableauLobattoIIIBIIIA(s)
    PartitionedTableau(Symbol("LobattoIIIBIIIA", s), TableauLobattoIIIB(s), TableauLobattoIIIA(s))
end

"Tableau for Gauss-Lobatto IIIC-IIIC̄ method with s stages"
function TableauLobattoIIICIIIC̄(s)
    PartitionedTableau(Symbol("LobattoIIICIIIC̄", s), TableauLobattoIIIC(s), TableauLobattoIIIC̄(s))
end

"Tableau for Gauss-Lobatto IIIC̄-IIIC method with s stages"
function TableauLobattoIIIC̄IIIC(s)
    PartitionedTableau(Symbol("LobattoIIIC̄IIIC", s), TableauLobattoIIIC̄(s), TableauLobattoIIIC(s))
end
