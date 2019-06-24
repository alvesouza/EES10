%Primeira questão%

%Definindo G e os requisitos do sistema
G = tf([500], [1 19 110 200 0])

Mp=.2;  %Overshoot
Tp=1.5; %Instante de pico 

%Calculando csi e wn
xi = -log(Mp)/sqrt(pi^2+[log(Mp)]^2)
wn  = pi/(Tp*sqrt(1-xi^2))

%Calculando a margem de fase e wc
PMd = 100*xi %Verifica-se que a margem de fase desejada é menor que a margem de fase na freuência wc -> projete um controlador de avanço, pois somente um ganho k não é suficiente
wc = wn

%Realiza bode para achar o a fase e o módulo em G(j*wc) 
[m, p] = bode(G,[wc])
%Ganho para atender o requisito de tempo de subida
k = 1/m;
kG = series(k, G); margin(kG) %Verifica-se que a margem de fase desejada não é atendida somente com um ganho K

%Calcula margem de fase atual
PMa = 180+p
%Calcula o avançoo, que é a margem de fase desejada menos a atual
Av = PMd-PMa %Av = PMd-180-p

%Calcula o alfa característico do avançoo de fase
alfa  = (1-sind(Av))/(1+sind(Av))
%Calcula a função de trasnferencia
Cfreq = zpk([-sqrt(alfa)*wc],[-wc/sqrt(alfa)],[1/m/sqrt(alfa)])

%Calcula a funÃ§Ã£o de transferÃªncia de G em serie com o controlador
CGfreq = series(Cfreq,G); margin(G*Cfreq) %Verifica-se que margem de fase foi atendida

%Simulando a resposta ao degrao, e verificando se atende aos requisitos
%necessários
stepinfo(feedback(CGfreq, 1))

%Verifica-se que o overshooot é maior que o desejado. Admitindo uma redução
%de 10% no ganho do controlador tem-se

Cfreq = zpk([-sqrt(alfa)*wc],[-wc/sqrt(alfa)],[0.9*1/m/sqrt(alfa)])
CGfreq = series(Cfreq,G)

%Simulando novamente:
stepinfo(feedback(CGfreq, 1))


%FunÃ§Ã£o de trasnferÃªncia ideal com os polos desejados
Td = tf([wn^2],[1 2*xi*wn wn^2])
%SoluÃ§Ã£o em resposta em frequÃªncia
Tfreq = feedback(CGfreq, 1)


%Step de todas as soluÃ§Ãµes
figure; step(Td, Tfreq); legend TOGGLE
