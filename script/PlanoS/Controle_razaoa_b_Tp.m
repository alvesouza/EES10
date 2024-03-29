function C = Controle_razaoa_b_Tp(G,razaoa_b,Mp,Tp,T1,T2,T3)
%Mp � o overshoot
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
fprintf("q = x+i*y\n");
fprintf("(K/abs(q+b))*(q+razaoa_b*b)*(x+b - j*y)\n\n");
fprintf("Encontrando a fase do Controlador\n");
fprintf("x^2 +(r+1)*x*b+r*b^2+y^2+(1-r)*y*b*j=");
fprintf("(x^2 +(r+1)*x*b+r*b^2+y^2)*tan(Phi) = (1-r)*y*b");
fprintf("razaoa_b*tand(Phi)*b^2 + ((razaoa_b+1)*tand(Phi)*real(q) +real(q)*(razaoa_b-1))*b + tand(Phi)*(abs(q)^2) = 0\n")
%Calculando o polo e a raiz do controlador
pol = [razaoa_b*tand(Phi) ((razaoa_b+1)*tand(Phi)*real(q) +(razaoa_b-1)*imag(q)) tand(Phi)*(abs(q)^2)];
fprintf("%f*b^2 + (%f)*b + (%f) = 0\n",pol(1),pol(2),pol(3));
b = roots(pol);
fprintf("a = razaoa_b*b\n");
a = razaoa_b*b;
K = K*[abs(q+b(1))/abs(q+a(1)) abs(q+b(2))/abs(q+a(2))];
fprintf("K = %s\n",mat2str(K));
fprintf("a = %s\n",mat2str(a));
fprintf("b = %s\n",mat2str(b));
C = [zpk([-a(1)],[-b(1)],[K(1)]) zpk([-a(2)],[-b(2)],[K(2)])];
end