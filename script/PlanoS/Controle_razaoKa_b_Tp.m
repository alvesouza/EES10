function C = Controle_razaoKa_b_Tp(G,razaoKa_b,Mp,Tp,T1,T2,T3)
%Mp é o overshoot
%razaoKa_b valor desejado de k*a/b do controlador
%Tp tempo de pico
%T1 retardo antes da malha fechada
%T2 retardo antes da malha malha direta
%T3 retardo antes da malha malha de retorno
[z,p,k] = zpkdata(G);
raizes = cell2mat(z);
polos = cell2mat(p);

fprintf("Determinando os valores de G\n");
fprintf("k = %d\n",k);
fprintf("raizes = %s\n",mat2str(raizes));
fprintf("polos = %s\n\n",mat2str(polos));

j = sqrt(-1);
nroot = size(raizes,1);
npoles = size(polos,1);
K = 1;
%Calculando csi e wn
xi = -log(Mp)/sqrt(pi^2+[log(Mp)]^2);
wn  = pi/((Tp-T1-T2)*sqrt(1-xi^2));

%Plano-s
q = -xi*wn + j*wn*sqrt(1-xi^2);
fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("xi = %f\nwn = %f\n",xi,wn);
fprintf("poloDesejado = %f + j*%f\n\n",real(q),imag(q));

Phi = phase(k);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("Phi é o valor da fase que o Controlador deve possuir\n");
fprintf("Phi = -180 + phased(k) - fase(raizes) + fase(polos) + imag(q)*(T2+T3)\n");

fprintf("Se k negativo ou feedback positivo,");
fprintf("a formula de Phi não possui -180\n\n");

fprintf("phased(k) = %f\n\n",Phi*180/pi);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("Calculando as fases das raizes\n");
fase = 0;
i = 1;
while i<= nroot
    
    fprintf("raiz = %f + i*(%f)\n",real(raizes(i)),imag(raizes(i)));
    fprintf("fase(q-raizes(i))=>fase(%f+j*(%f)) = %f\n",real(q-raizes(i)),imag(q-raizes(i)),phase(q-raizes(i))*180/pi);
    
    Phi = Phi-phase(q-raizes(i));
    
    fase = fase-phase(q-raizes(i));
    i = i + 1;
end
fprintf("\nfase(q-raizes) = %f\n\n",fase*180/pi);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("Calculando as fases dos polos\n");

fase = 0;
i = 1;
while i<= npoles
    
    fprintf("polo = %f + i*%f\n",real(polos(i)),imag(polos(i)));
    fprintf("fase(q-polos(i))=>fase(%f+j*(%f)) = %f\n",real(q-polos(i)),imag(q-polos(i)),phase(q-polos(i))*180/pi);
    
    Phi = Phi+phase(q-polos(i));
    
    fase = fase+phase(q-polos(i));
    
    i = i + 1;
end

fprintf("\nfase(q-polos) = %f\n",fase*180/pi);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("\nimag(q)*(T2+T3)*180/pi = %f*(%f+%f)*180/pi => %f\n",imag(q),T2,T3,imag(q)*(T2+T3)*180/pi);
Phi = Phi+imag(q)*(T2+T3);
Phi = Phi/pi*180;
Phi = -180  + Phi;
Phi = phase(cosd(Phi) + j*sind(Phi))*180/pi;

fprintf("\nO valor da fase do Controlador em q\n");
fprintf("\nPhi = %f\n\n\n",Phi);

%calculando o valor do K

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("Calculando o valor de K,que o controlador deve ter.\n\n");
fprintf("K*abs(q+a)/abs(q+b) = (abs(q-polos))/(abs(k)*abs(q-raizes)*");
fprintf("exp(-real(q)*(T2+T3)))\n\n");

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("Calculando o valor de K*abs(q+a)/abs(q+b)\n\n");

fprintf("Calculando os abs(q-polos)\n");
fase = 1;
i = 1;
while i<= npoles
    
    fprintf("polo = %f + i*%f\n",real(polos(i)),imag(polos(i)));
    fprintf("abs(q-polos(i))=>abs(%f + j*(%f)) = %f\n",real(q-polos(i)),imag(q-polos(i)),abs(q-polos(i)));
        
    K = K*abs(q-polos(i));
    fase = fase*abs(q-polos(i));
    i = i + 1;
end

fprintf("\nabs(q-polos) = %f\n\n\n",fase);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("\nCalculando os abs(q-raizes)\n");
fase = 1;
i = 1;
while i<= nroot
    
    fprintf("raiz = %f + i*%f\n",real(raizes(i)),imag(raizes(i)));
    fprintf("abs(q-raizes(i))=>abs(%f + j*(%f)) = %f\n",real(q-raizes(i)),imag(q-raizes(i)),abs(q-raizes(i)));

    K = K/abs(q-raizes(i));
    fase = fase*abs(q-raizes(i));
    i = i + 1;
end

fprintf("\nabs(q-raizes) = %f\n\n",fase);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("abs(k)*exp(-real(q)*(T2+T3))=>abs(%f)*exp(-%f*(%f+%f))= %f\n",k,real(q),T2,T3,abs(k)*exp(-real(q)*(T2+T3)));


K = K/(abs(k)*exp(-real(q)*(T2+T3)));
fprintf("\nK*abs(q+a)/abs(q+b) = %f\n\n\n",K);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("Encontrando os valores do Controlador\n");
fprintf("Fazendo beta = 1/|G(q)| = %f;\n",K);
beta = K;
fprintf("q = x+i*y\n");
fprintf("Temos a relação:\n\n");
fprintf("k(q+a)/(q+b) =|k(q+a)/(q+b)|(cos(Phi)+j*sin(Phi)) =");
fprintf("beta(cos(Phi)+j*sin(Phi))\n=>");
fprintf("q*k+(razaoKa_b-beta*(cos(Phi)+j*sin(Phi)))*b = q*beta*(cos(Phi)+j*sin(Phi)))\n=>");
fprintf("Substituindo os valores, para encontrar as matrizes para resolver a equação de 2 variaveis")
fprintf("%s*k+(%f-%f*(%f+j*(%f)))*b = %s*%f*(%f+j*(%f)))\n=",q,razaoKa_b,beta,cosd(Phi),sind(Phi),beta,cosd(Phi),sind(Phi));

vetorEsq = [q razaoKa_b-beta*(cosd(Phi)+j*sind(Phi))]%corresponde a k e b
iqualdadeDir = q*beta*(cosd(Phi)+j*sind(Phi))%valor da igualdade da equação

fprintf("Como os valores de k e b são reais, logo podemos dividir a equação acima em real e imaginario:\n\n")
realVetorEsq = real(vetorEsq)
imagVetorEsq = imag(vetorEsq)
realIqualdadeDir = real(iqualdadeDir)
imagIqualdadeDir = imag(iqualdadeDir)

fprintf("temos AX = B=>\n")
A = [realVetorEsq(1) realVetorEsq(2);imagVetorEsq(1) imagVetorEsq(2)]
fprintf("Para A primeira linha corresponde a parte real da equação esquerda e coluna ao multiplicativo de k")
B = [realIqualdadeDir;imagIqualdadeDir]
fprintf("Para B primeira linha corresponde a parte real da equação direita")
fprintf("X = inv(A)*B=>\n")
X = A\B%inv(A)*B só que mais rapido e preciso
k = X(1);
b = X(2);
fprintf("k = %f\nb = %f\n",k,b);
a = razaoKa_b*b/k;
C = zpk([-a],[-b],[k]);
end