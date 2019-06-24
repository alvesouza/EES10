function k = k_Gain(G,Mg)
[m,p] = margin(G);
kf = 10^(Mg/20);
k = m/kf;
end