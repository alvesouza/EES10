j = sqrt(-1);
T1 = 0;
T2 = 0;
T3 = 0;
erro = 0.05;
s = tf('s');
delay1 = exp(-T1*s);%fora
delay2 = exp(-T2*s);%malha direta
delay3 = exp(-T3*s);%feedback

%dividindo os delay equivalentes
delayout = pade(delay1*delay2,1);
delayfeed = pade(delay2*delay3,1);
Mp = 0.095;
xi = -log(Mp)/sqrt(pi^2+[log(Mp)]^2);
%%
%Tempo de pico
Tp = 0.8;
wn  = (pi)/((Tp-T1-T2)*sqrt(1-xi^2));
%Polo
q = -xi*wn + j*wn*sqrt(1-xi^2);
%%
%Tempo de subida
Tr = 3;
wn  = (pi-acos(xi))/((Tr-T1-T2)*sqrt(1-xi^2));
%Polo
q = -xi*wn + j*wn*sqrt(1-xi^2);
%%
%fazendo por tf
dem = [1 19 110 200 0];
num = [500];
k = num(1)/dem(1);
polos = roots(dem)';
raizes = roots(num)';
G = tf(num, dem);
%G = zpk(raizes,polos,k)

C = Controle_equi_Tp(k,raizes,polos,Mp,Tp,T1,T2,T3);

%separando as raizes,polos e proporcional do Controle
[z,p,k] = zpkdata(C);
z = cell2mat(z);
p = cell2mat(p);

%modificando os valores encontrados para satisfazer as condições
C = zpk(1.0*z,1.5*p,2.0*k);
T = delayout*feedback(C*G,delayfeed);
%copie e cole no terminal o abaixo para encontrar o apropriado
stepinfo(delayout*feedback(zpk(1.0*z,0.98*p,1.125*k)*G,delayfeed))
%step(T)

%%
%fazendo por zpk
K = 10;
polos = [0 -5 -10];
raizes = [5];
G = zpk(raizes,polos,K);
C = Controle_a_Tr(K,raizes,polos,0.25,Mp,Tr,T1,T2,T3)

%separando as raizes,polos e proporcional do Controle
[z,p,k] = zpkdata(C);
z = cell2mat(z);
p = cell2mat(p);

%modificando os valores encontrados para satisfazer as condições
C = zpk(1.0*z,1.0*p,1.0*k);
T = delayout*feedback(C*G,delayfeed);
%copie e cole no terminal o abaixo para encontrar o apropriado
stepinfo(delayout*feedback(zpk(1.0*z,1.0*p,1.0*k)*G,delayfeed),'RiseTimeLimits',[0,1])
%stepinfo(feedback(G*Controle_a_Tr(K,raizes,polos,0.25,Mp,Tr,T1,T2,T3),1),'RiseTimeLimits',[0,1])
stepinfo(feedback(Controle_b_Tp(K,raizes,polos,-5,Mp,Tp,T1,T2,T3)*G3,1))