F = zpk([-5],[-1],[1])
G = zpk([],[-3],[1])
H = zpk([],[0],[1])

T1 = series(G,H)
T2 = feedback(T1,1,+1)
T = parallel(T2,F)

step(T)

rlocus(-T1)