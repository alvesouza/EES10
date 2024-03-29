function C = Controle_equi_Tr(G,Mp,Tr,T1,T2,T3)
%Mp � o overshoot
%Tr tempo de subida 0-100%
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
wn  = (pi-acos(xi))/((Tr-T1-T2)*sqrt(1-xi^2));

%Plano-s
q = -xi*wn + j*wn*sqrt(1-xi^2);
fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("xi = %f\nwn = %f\n",xi,wn);
fprintf("poloDesejado = %f + j*%f\n\n",real(q),imag(q));

Phi = phase(k);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("Phi � o valor da fase que o Controlador deve possuir\n");
fprintf("Phi = -180 + phased(k) - fase(raizes) + fase(polos) + imag(q)*(T2+T3)\n");

fprintf("Se k negativo ou feedback positivo,");
fprintf("a formula de Phi n�o possui -180\n\n");

fprintf("phased(k) = %f\n\n",Phi*180/pi);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("Calculando as fases das raizes\n");
fase = 0;
i = 1;
while i<= nroot
    
    fprintf("raiz = %f + i*(%f)\n",real(raizes(i)),imag(raizes(i)));
    fprintf("fase(q-raizes(i))=>fase(%f+j*(%f)) = %f\n",real(q-raizes(i)),imag(q-raizes(i)),phase(q-raizes(i))*180/pi);
    
    Phi = Phi-phase(q-raizes(i));
    
    fase = fase+phase(q-raizes(i));
    i = i + 1;
end
fprintf("\nfase(q-raizes) = %f\n\n",fase*180/pi);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||");
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

%Calculando o polo e a raiz do controlador
fprintf("fase(q+a) = 90 + Phi/2; fase(q+b) = 90-Phi/2\n");
fprintf("Phi/2 = %f\n",Phi/2);
fprintf("a = -real(q) + imag(q)/tand(Phi/2)=%f + (%f)/tan(%f)\n",-real(q),imag(q),90+Phi/2);
fprintf("b = -real(q) + imag(q)/tand(-Phi/2)=%f + (%f)/tan(%f)\n",-real(q),imag(q),90-Phi/2);
a = -real(q) + imag(q)/tand(90 + Phi/2);
b = -real(q) + imag(q)/tand(90 - Phi/2);
fprintf("a = %f\n",a);
fprintf("b = %f\n",b);

fprintf("K = %f*abs(%f + i*(%f))/abs(%f + i*(%f))\n\n",K,real(q+b),imag(q+b),real(q+a),imag(q+a));
fprintf("K = %f\n",K);
fprintf("a = %f\n",a);
fprintf("b = %f\n",b);

C = zpk([-a],[-b],[K]);
end