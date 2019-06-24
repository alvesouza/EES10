function C = Lead_LagS_a_Tr(G,a,Mp,e,Tr,T1,T2,T3)
[z,p,k] = zpkdata(G);
raizes = cell2mat(z);
polos = cell2mat(p);
Kpo = k;
j = sqrt(-1);
razao = 100;
%Calculando csi e wn
xi = -log(Mp)/sqrt(pi^2+[log(Mp)]^2);
wn  = (pi-acos(xi))/((Tr-T1-T2)*sqrt(1-xi^2));

%Plano-s
q = -xi*wn + j*wn*sqrt(1-xi^2);
Cav = Controle_equi_Tr(G,a,Mp,Tr,T1,T2,T3);

nroot = size(raizes,1);
npoles = size(polos,1);
[z,p,k] = zpkdata(Cav);
z = cell2mat(z);
p = cell2mat(p);
fprintf("\nTesmos o Kpf = 1/e - 1 = 1/%f - 1\n",e);
Kpf = 1/e - 1;
fprintf("Kpf = %f\n",Kpf);
Kpo = Kpo*k;

i = 1;
while i <= nroot
    Kpo = Kpo*(-raizes(i));
    i = i + 1;
end
i = 1;
while i <= npoles
    Kpo = Kpo/(-polos(i));
    i = i + 1;
end
Kpo = Kpo*(z/p);
a_b = Kpf/Kpo;
fprintf("a/b = %f\n",a_b);
if abs(a_b) > 1
    a = real(q)/razao;
    b = a/a_b;
else
    b = real(q)/razao;
    a = b*a_b;
end
if isnan(a)
    a = 0;
end
if isnan(b)
   b = 0; 
end
C = Cav*zpk([a],[b],[1]);
end