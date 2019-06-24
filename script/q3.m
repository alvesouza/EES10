G1 = zpk([],[0 0],[1])
G2 = zpk([],[-1],[50])
G3 = zpk([],[0],[2])
G4 = zpk([0],[],[1])
G5 = zpk([],[],[-2])

Aux1 = feedback(G2,G3)
Aux2 = parallel(G4,-G5)
Aux3 = series(Aux1,Aux2)
Aux4 =series(G1,Aux3)

T = feedback(Aux4,1)

step(T)
