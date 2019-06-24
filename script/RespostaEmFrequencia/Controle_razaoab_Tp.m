function C = Controle_razaoab_Tp(G,razaoab,PMd,Tp)

[z,p,k] = zpkdata(G);
raizes = cell2mat(z);
polos = cell2mat(p);

fprintf("Determinando os valores de G\n");
fprintf("k = %d\n",k);
fprintf("raizes = %s\n",mat2str(raizes));
fprintf("polos = %s\n\n",mat2str(polos));

%caso seja de atraso
[m,p] = margin(G);
PMa = 180+p;
fprintf("PMa = %f\n",PMa);
if PMa>PMd
    PMd = PMd + 5;
    fprintf("Logo novo PMd = %f\n",PMd);
end
j = sqrt(-1);
%Calculando csi e wn

xi = PMd/100;
wc  = pi/(Tp*sqrt(1-xi^2));
fprintf("Calculando o wc = %f\n",wc);
[m,p] = margin(G,[wc]);
PMa = 180+p-phase(k)*180/pi;
fprintf("novo PMa = %f\n",PMa);
fprintf("Tem Av = Pmd-Pma = %f - %f\n",PMd,PMa);
Av = PMd-PMa; %Av = PMd-180-p
fprintf("Av = %f\n",Av);

denominador = conv([1 j*wc],[1 -j*wc]);
numerador = conv([razaoab j*wc],[1 -j*wc]);
fprintf("G(jwc) = (K/(dem*dem'))*(num*dem') => (K/(dem*dem'))*A\n")
fprintf("Como Polinomio de b:\nA = %s + j*%s",mat2str(real(numerador)),mat2str(imag(numerador)));
fprintf("faseA = Av = %f\n",Av);
fprintf("ImA = RealA*tang(Av)\n")
pol = tand(Av)*real(numerador) - imag(numerador);
fprintf("Polinomio de b\n");
fprintf("%s = 0\n",mat2str(pol));
b = roots(pol);
a = razaoab*b;
fprintf("b = %s\n",mat2str(b));
fprintf("a = %s\n",mat2str(a));
fprintf("Para encontrar o k\nk=1/m*abs(b+wc*j)/abs(a+wc*j)")
k = [1/m*abs(b(1)+wc*j)/abs(a(1)+wc*j) 1/m*abs(b(2)+wc*j)/abs(a(2)+wc*j)];
fprintf("k = %s",mat2str(k));
fprintf("a = %s\n",mat2str(a));
fprintf("b = %s\n",mat2str(b));
C=[zpk([-a(1)],[-b(1)],[k(1)]) zpk([-a(2)],[-b(2)],[k(2)])];

end