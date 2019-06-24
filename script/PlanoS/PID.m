function C = PID(G,Kv,Mp,Tp,T1,T2,T3)
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
x = real(q);
y = imag(q);
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
    
    fase = fase+phase(q-raizes(i));
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
fprintf("Calculando o valor de |G(q)|\n\n");

fprintf("Calculando os abs(q-polos)\n");
fase = 1;
i = 1;
while i<= npoles
    
    fprintf("polo = %f + i*%f\n",real(polos(i)),imag(polos(i)));
    fprintf("abs(q-polos(i))=>abs(%f + j*(%f)) = %f\n",real(q-polos(i)),imag(q-polos(i)),abs(q-polos(i)));
        
    K = K/abs(q-polos(i));
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

    K = K*abs(q-raizes(i));
    fase = fase*abs(q-raizes(i));
    i = i + 1;
end

fprintf("\nabs(q-raizes) = %f\n\n",fase);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("abs(k)*exp(-real(q)*(T2+T3))=>abs(%f)*exp(-%f*(%f+%f))= %f\n",k,real(q),T2,T3,abs(k)*exp(-real(q)*(T2+T3)));


K = K*(abs(k)*exp(-real(q)*(T2+T3)));
fprintf("\n|G(q)| = %f\n\n\n",K);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");

%Calculando o polo e a raiz do controlador
fprintf("Calculando a limG\n");
limG = Gval(G,T1,T2,T3,0);
fprintf("Calculando a Ki = Kv/limG\n");
Ki = Kv/limG;
fprintf("Ki = %f\n\n",Ki);

%Calculando a equação da relação da fase
fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("alpha = <(Kd*q^2+Kp*q+Ki)= -180-phi+<(q) = -180-(%f)+(%f)",Phi,phase(q)*180/pi);
alpha = -180-Phi+phase(q)*180/pi;
fprintf("= %f\n\n",alpha);
Eq = [abs(q)^2-2*x*y/tand(alpha) x+y/tand(alpha);0 0];
Resul = [-Ki -Ki];
fprintf("(abs(q)^2-2*x*y/tand(alpha))*Kd + (x+y/tand(alpha))*Kp = -Ki\n");
fprintf("(%f)*Kd + (%f)*Kp = %f\n",Eq(1,1),Eq(1,2),Resul(1,1));
fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
A = abs(q)*cosd(alpha)/K;
fprintf("|Kd*abs(q)^2+Kp*x+Ki +i(Kp*y+2*x*y*Kd)|=|q|/|G(q)|\n");
fprintf("Kd*abs(q)^2+Kp*x+Ki = +-|q|*cosd(alpha)/|G(q)|\n");
Eq(2,1) = abs(q)^2;
Eq(2,2) = x;
fprintf("Kd*(%f)^2+Kp*(%f) = %f +-(%f)\n",Eq(2,1),Eq(2,2),-Ki,A);
fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
Resul1 = Resul;
Resul2 = Resul;
Resul1(2) = Resul1(2) + A;
Resul2(2) = Resul2(2) - A;
fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("Resolvendo as equaçãoes\n");
Sol1 = (Eq\Resul1')';
Sol2 = (Eq\Resul2')';
fprintf("Solução(+) 1\n\n");
fprintf("Kd = %f;Kp = %f;Ki = %f\n",Sol1(1),Sol1(2),Ki);
fprintf("Solução(-) 2\n\n");
fprintf("Kd = %f;Kp = %f;Ki = %f\n",Sol2(1),Sol2(2),Ki);
C = [tf([Sol1(1) Sol1(2) Ki],[1 0]) tf([Sol2(1) Sol2(2) Ki],[1 0])]
end