%%Q1
G1 = zpk([],[-1],[2])
G2 = zpk([-2],[-3],[1])
G3 = zpk([],[0],[2])

Aux1 = feedback(G1,G2)
Aux2 = series(Aux1,G3)

T = feedback(Aux2,1)
damp(T)
stepinfo(T)
step(T)

%%
%%Q2
T1 =feedback(-G1, G2)
T2 =feedback(G1, -G2)
%%
%%Q3

G1 = zpk([-1],[-2 -5 -10],[100])
G2 = tf([64], [1 6 64])