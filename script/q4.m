G1 = zpk([],[-5],[10])
G2 = zpk([],[0],[1])
G3 = zpk([],[],[5])

Aux1 = feedback(G1,G3)

Aux2 = series(Aux1,G2)

T = feedback(Aux2,1)

step(T)

rlocus(Aux2)

bode(Aux2)

%margem de fase 90 graus
%margem de ganho inf
