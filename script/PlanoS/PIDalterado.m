function C = PIDalterado(G,Kv,Mp,Tr,T1,T2,T3)
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
x = real(q);
y = imag(q);
fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("xi = %f\nwn = %f\n",xi,wn);
fprintf("poloDesejado = x+y*j = %f + j*%f\n\n",real(q),imag(q));

Phi = -180;

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("Phi é o valor da fase que o Controlador deve possuir\n");
fprintf("Phi = -180 - phased(G(q))\n");


fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("Calculando o valor de G(q)\n");
Gq = Gval(G,T1,T2,T3,q);
fprintf("G(q) = %f + (%f)*j\n\n",real(Gq),imag(Gq));
fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("Calculando o valor de |G(q)|\n\n");
modG = abs(Gq);
fprintf("|G(q)| = %f\n",modG);
fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("Calculando o valor de phased(G(q))\n");
phasedGq = phase(Gq)*180/pi;
fprintf("phased(G(q)) = %f\n",phasedGq);
Phi = -180 - phasedGq;
fprintf("Phi = %f\n\n",Phi);
fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
%Calculando o polo e a raiz do controlador
fprintf("Calculando a limG\n");
limG = Gval(G,T1,T2,T3,0);
fprintf("Calculando a Ki = Kv/limG\n");
Ki = Kv/limG;
fprintf("Ki = %f\n\n",Ki);

%Calculando a equação da relação da fase
fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("alpha = <(Kd*q^2+Kp*q+Ki)= phi+<(q) = (%f)+(%f)",Phi,phase(q)*180/pi);
alpha = Phi+phase(q)*180/pi;
fprintf("= %f\n\n",alpha);
Eq = [abs(q)^2-2*x*y/tand(alpha) x-y/tand(alpha);0 0];
Resul = [-Ki -Ki];
fprintf("(abs(q)^2-2*x*y/tand(alpha))*Kd + (x-y/tand(alpha))*Kp = -Ki\n");
fprintf("(%f)*Kd + (%f)*Kp = %f\n",Eq(1,1),Eq(1,2),Resul(1,1));
fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
A = abs(q)*cosd(alpha)/modG;
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