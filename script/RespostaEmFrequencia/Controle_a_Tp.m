function C = Controle_a_Tp(G,a,PMd,Tp)


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

PMa = 180+p;
fprintf("novo PMa = %f\n",PMa);
fprintf("Tem Av = Pmd-Pma = %f - %f\n",PMd,PMa);
Av = PMd-PMa; %Av = PMd-180-p
fprintf("Av = %f\n",Av);

denominador = conv([1 j*wc],[1 -j*wc]);
numerador = conv([1 -j*wc],[a+j*wc]);
fprintf("G(jwc) = (K/(dem*dem'))*(num*dem') => (K/(dem*dem'))*A\n")
fprintf("Como Polinomio de b:\nA = %s + j*%s",mat2str(real(numerador)),mat2str(imag(numerador)));
fprintf("faseA = Av = %f\n",Av);
fprintf("ImA = RealA*tang(Av)\n")
pol = tand(Av)*real(numerador) - imag(numerador);
fprintf("Polinomio de b\n");
fprintf("%s = 0\n",mat2str(pol));
b = roots(pol);
fprintf("b = %s\n",mat2str(b));
fprintf("Para encontrar o k\nk=1/m*abs(b+wc*j)/abs(a+wc*j)")
k = 1/m*abs(b+wc*j)/abs(a+wc*j);
fprintf("k = %f",k);
fprintf("a = %f\n",a);
fprintf("b = %s\n",mat2str(b));
C=zpk([-a],[-b(1)],[k]);

end