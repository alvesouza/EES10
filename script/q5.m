G1 = zpk([],[1 -9],[9])

figure
bode(G1)
figure
nyquist(G1)
figure
nichols(G1)
figure
rlocus(G1)


