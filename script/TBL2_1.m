%Primeira quest�o%

%Definindo G e os requisitos do sistema
G = tf([500], [1 19 110 200 0])

Mp=.2;  %Overshoot
Tp=1.5; %Instante de pico 

%Calculando csi e wn
xi = -log(Mp)/sqrt(pi^2+[log(Mp)]^2)
wn  = pi/(Tp*sqrt(1-xi^2))

%Calculando a margem de fase e wc
PMd = 100*xi %Verifica-se que a margem de fase desejada � menor que a margem de fase na freu�ncia wc -> projete um controlador de avan�o, pois somente um ganho k n�o � suficiente
wc = wn

%Realiza bode para achar o a fase e o m�dulo em G(j*wc) 
[m, p] = bode(G,[wc])
%Ganho para atender o requisito de tempo de subida
k = 1/m;
kG = series(k, G); margin(kG) %Verifica-se que a margem de fase desejada n�o � atendida somente com um ganho K

%Calcula margem de fase atual
PMa = 180+p
%Calcula o avan�oo, que � a margem de fase desejada menos a atual
Av = PMd-PMa %Av = PMd-180-p

%Calcula o alfa caracter�stico do avan�oo de fase
alfa  = (1-sind(Av))/(1+sind(Av))
%Calcula a fun��o de trasnferencia
Cfreq = zpk([-sqrt(alfa)*wc],[-wc/sqrt(alfa)],[1/m/sqrt(alfa)])

%Calcula a função de transferência de G em serie com o controlador
CGfreq = series(Cfreq,G); margin(G*Cfreq) %Verifica-se que margem de fase foi atendida

%Simulando a resposta ao degrao, e verificando se atende aos requisitos
%necess�rios
stepinfo(feedback(CGfreq, 1))

%Verifica-se que o overshooot � maior que o desejado. Admitindo uma redu��o
%de 10% no ganho do controlador tem-se

Cfreq = zpk([-sqrt(alfa)*wc],[-wc/sqrt(alfa)],[0.9*1/m/sqrt(alfa)])
CGfreq = series(Cfreq,G)

%Simulando novamente:
stepinfo(feedback(CGfreq, 1))


%Função de trasnferência ideal com os polos desejados
Td = tf([wn^2],[1 2*xi*wn wn^2])
%Solução em resposta em frequência
Tfreq = feedback(CGfreq, 1)


%Step de todas as soluções
figure; step(Td, Tfreq); legend TOGGLE
