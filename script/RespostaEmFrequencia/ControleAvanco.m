function C = ControleAvanco(G,PMd,Tp)

[z,p,k] = zpkdata(G);
raizes = cell2mat(z);
polos = cell2mat(p);

fprintf("Determinando os valores de G\n");
fprintf("k = %d\n",k);
fprintf("raizes = %s\n",mat2str(raizes));
fprintf("polos = %s\n\n",mat2str(polos));

j = sqrt(-1);
%Calculando csi e wn
xi = PMd/100;
wc  = pi/(Tp*sqrt(1-xi^2));
fprintf("Calculando o wc = %f\n",wc);
%[m,p] = margin(G,wc);
valG = evalfr(G,wc*j);
m = abs(valG);
p = phase(valG)*180/pi;
fprintf("m = |G(wc*j)| = %f\np = <G(wc*j)=%f",m,p);
fprintf("PMa = 180+<G(wc*j)-phase(k)= 180+(%f)-phase(%f) = ",p,k);
PMa = 180+p-phase(k)*180/pi;
fprintf("PMa = %f",PMa);
fprintf("Tem Av = Pmd-Pma = %f - %f\n",PMd,PMa);
Av = PMd-PMa; %Av = PMd-180-p
fprintf("Av = %f\n",Av);

alfa = (1-sind(Av))/(1+sind(Av));

fprintf("alfa = (1-sind(Av))/(1+sind(Av))\n")
fprintf("xi = %f;alfa = %f\n",xi,alfa)
fprintf("Cav = (1/(m*sqrt(alfa)))*(s+wc*sqrt(alfa))/(s+wc/sqrt(alfa))\n")
C=zpk([-sqrt(alfa)*wc],[-wc/sqrt(alfa)],[1/m/sqrt(alfa)]);

end