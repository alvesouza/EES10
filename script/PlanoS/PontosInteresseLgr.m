function [alpha,angulos,raizesDerivada] = PontosInteresseLgr(G)

j = sqrt(-1);

%separando a função G
[z,p,k] = zpkdata(G);
raizes = cell2mat(z);
polos = cell2mat(p);
fprintf("Determinando os valores de G\n");
fprintf("k = %d\n",k);
fprintf("raizes = %s\n",mat2str(raizes));
fprintf("polos = %s\n\n",mat2str(polos));
nroot = size(raizes,1);
npoles = size(polos,1);

alpha = 0;%encontro das abscissas.

i = 1;
while i<= nroot
    alpha = alpha + raizes(i);
    i = i + 1;
end

i = 1;
while i<= npoles
    alpha = alpha - polos(i);
    i = i + 1;
end

alpha = -alpha/abs(nroot-npoles);
fprintf("temos o ponto de encontro das assintotas\nalpha =(polos-zeros)/(npolos - nzeros)= %f\n",alpha)
angulo = -180*phase(-k)/(pi*abs(nroot-npoles));
dAngle = 360/abs(nroot-npoles);

angulos = angulo*1.^(0:(abs(nroot-npoles) - 1)) + dAngle*(0:(abs(nroot-npoles) - 1));
fprintf("Temos os angulos das assintotas\nangle = -180*phase(-k)/pi + 360/(nroot-npoles)\n");
fprintf("angulos = %s\n",mat2str(angulos));
%derivada de 1/K

%encontrando o numerador e a derivada do numerador 
numerador = [1];
i = 1;
while i<= nroot
    numerador = conv(numerador,[1 raizes(i)]);
    i = i + 1;
end

numeradorDer = polyder(numerador);

%encontrando o denominador e sua derivada
denominador = [1];
i = 1;
while i<= npoles
    denominador = conv(denominador,[1 polos(i)]);
    i = i + 1;
end

denominadorDer = polyder(denominador);
%[q,d] = polyder(a,b);
pol = conv(denominador,numeradorDer) - conv(denominadorDer,numerador);

raizesDerivada = roots(pol);
fprintf("Encontrando o ponto de encontro do ramo imaginario,com o real:\n");
fprintf("d(1/k) = d(1/G)/ds = 0\n");
fprintf("vetor do numerador:\n%s\n\n",mat2str(pol));
fprintf("as raizes da equação são: \n%s\n",mat2str(raizesDerivada))
%Routh-Hurwitz
fprintf("Faremos o Routh-Hurwitz na mão\n")
end