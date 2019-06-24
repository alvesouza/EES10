function C = Lead_LagS_equi_Kv_Tp(G,Mp,Kv,Tp,T1,T2,T3)
j = sqrt(-1);
razao = 10;

[z,p,k] = zpkdata(G);
raizes = cell2mat(z);
polos = cell2mat(p);
Kpo = k;

%Calculando csi e wn
xi = -log(Mp)/sqrt(pi^2+[log(Mp)]^2);
wn  = pi/((Tp-T1-T2)*sqrt(1-xi^2));

%Plano-s
q = -xi*wn + j*wn*sqrt(1-xi^2);
fprintf("Resolvendo o controlador primeiro\n");
Cav = Controle_equi_Tp(G,Mp,Tp,T1,T2,T3);

nroot = size(raizes,1);
npoles = size(polos,1);
[z,p,k] = zpkdata(Cav);
z = cell2mat(z);
p = cell2mat(p);

Kpf = Kv;
Kpo = Kpo*k;

i = 1;
while i <= nroot
    Kpo = Kpo*(-raizes(i));
    i = i + 1;
end
i = 1;
a = 1;
while i <= npoles
    if abs(polos(i)) > 0 || a==0
        Kpo = Kpo/(-polos(i));
    elseif a==1
    	a = 0;  
    end
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