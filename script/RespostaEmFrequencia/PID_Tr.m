function C = PID_Tr(G,PMd,Kv,Tr)

[z,p,k] = zpkdata(G);
raizes = cell2mat(z);
polos = cell2mat(p);
limG = k;
modG = k;

nroot = size(raizes,1);
npoles = size(polos,1);

i = 1;
while i <= nroot
    limG = limG*(-raizes(i));
    i = i + 1;
end
i = 1;
while i <= npoles
    limG = limG/(-polos(i));
    i = i + 1;
end

fprintf("Determinando os valores de G\n");
fprintf("k = %d\n",k);
fprintf("raizes = %s\n",mat2str(raizes));
fprintf("polos = %s\n\n",mat2str(polos));
fprintf("limG = %f\n",limG);
j = sqrt(-1);
%Calculando csi e wn

xi = PMd/100;
wc  = (pi-acos(xi))/(Tr*sqrt(1-xi^2));
fprintf("Calculando o wc = %f\n",wc);
fprintf("Calculando modG\n");

i = 1;
while i <= nroot
    modG = modG*abs(j*wc-raizes(i));
    i = i + 1;
end
i = 1;
while i <= npoles
    modG = modG/abs(j*wc-polos(i));
    i = i + 1;
end
fprintf("modG=%f\n",modG);
[m, p] = bode(G,[wc]);
PMa = 180+p;
Ki = Kv/limG;
fprintf("Ki=Kv/limG = %f/%f = %f\n",Kv,limG,Ki);
fprintf("PMa = %f\n",PMa);
fprintf("Tem Av = Pmd-Pma = %f - %f\n",PMd,PMa);
Av = PMd-PMa; %Av = PMd-180-p
fprintf("Av = %f\n",Av);
Kp = cosd(Av)/modG;

fprintf("Kp=cosd(Av)/modG = cosd(%f)/%f = %f\n",Av,modG,Kp);

Kd = sind(Av)/(modG*wc) + Ki/(wc^2);
fprintf("Kd = sind(Av)/(modG*wc) + Kt/(wc^2) = %f\n",Kd);

C=tf([Kd Kp Ki],[1 0]);

end