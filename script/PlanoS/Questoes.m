C1 = Controle_b_Tp(K,raizes,polos,4,Mp,Tp,T1,T2,T3)
P1 = zpk([-10],[-0.78],[1])

C3 = Controle_b_Tp(K,raizes,polos,3,Mp,Tp,T1,T2,T3)
P1 = zpk([-10],[-0.78],[1])

P = tf([1 12.31],[1 -5])
ganhoFeed = -1.2;