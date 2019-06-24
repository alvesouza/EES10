G1 = zpk([],[-2],[2])
G2 = zpk([],[0],[1])
G3 = zpk([],[],[5])
G4 = zpk([-5],[-3, -10],[1])
G5 = zpk([],[],[2])
G6 = zpk([0],[],[1])
G7 = zpk([],[0],[5])

Aux1 = feedback(G1,G2)
Aux2 = feedback(G4,G5)
Aux3 = series(Aux1,G3)
Aux4 = series(Aux2,G6)
Aux5 = series(Aux3,Aux4)
Aux6 = series(Aux3,Aux2)

T = parallel(Aux5,G7)

step(Aux6)