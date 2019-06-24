%%close all
%clear all
%clc
G=zpk([],[0 -2 -10],[50])
Mp=.15; 
Tr=1.2; 

xi=-log(Mp)/sqrt(pi^2+[log(Mp)]^2)
wn  = (pi-acos(xi))/Tr/sqrt(1-xi^2)

PMd = 100*xi
wc = wn

[m, p] = bode(G,[wc])
k = 1/m;
kG = series(k, G); margin(kG)

PMa = 180+p
Av = PMd-PMa %Av = PMd-180-p

alfa = (1-sind(Av))/(1+sind(Av))
Cfreq=zpk([-sqrt(alfa)*wc],[-wc/sqrt(alfa)],[1/m/sqrt(alfa)])

CGfreq = series(Cfreq,G); margin(G*Cfreq)

bode(kG, CGfreq); legend TOGGLE
figure; nichols(G, kG, CGfreq); grid on;  legend TOGGLE

%Plano-s
q = -xi*wn + j*wn*sqrt(1-xi^2)
PGq = (-phase(q+2)-phase(q+10)-phase(q))/pi*180

Phi = -PGq-180

ar = -real(q)
br = -real(q) + imag(q)/tand(270+PGq)

aer = -real(q) -imag(q)*tand((-PGq-180)/2)
ber = -real(q) +imag(q)*tand((-PGq-180)/2)

(phase(q+ar)-phase(q+br))/pi*180
(phase(q+aer)-phase(q+ber))/pi*180


kr = (abs(q+2)*abs(q+10)*abs(q)*abs(q+br))/(50*abs(q+ar))
ker = (abs(q+2)*abs(q+10)*abs(q))/(50)

Cr = zpk([-ar],[-br],[kr])
Cer = zpk([-aer],[-ber],[ker])

CrG = series(Cr, G)
CerG = series(Cer, G)

rlocus(CrG, CerG); legend TOGGLE

Td=tf([wn^2],[1 2*xi*wn wn^2])
Tfreq=feedback(CGfreq, 1)
Tr=feedback(CrG, 1)
Ter=feedback(CerG, 1)

damp(Td)
damp(Tr)
damp(Ter)

step(Td, Tfreq, Tr, Ter); legend TOGGLE


stepinfo(Td,'RiseTimeLimits',[0,1])
stepinfo(Tfreq,'RiseTimeLimits',[0,1])
stepinfo(Tr,'RiseTimeLimits',[0,1])
stepinfo(Ter,'RiseTimeLimits',[0,1])

Tfreq1=feedback(.9*CGfreq, 1)
Tr1=feedback(1.1*CrG, 1)
Ter1=feedback(1.1*CerG, 1)
step(Td, Tfreq1, Tr1, Ter1); legend TOGGLE

Tfreq2=feedback(.8*CGfreq, 1)
Tr2=feedback(1.2*CrG, 1)
Ter2=feedback(1.2*CerG, 1)
step(Td, Tfreq2, Tr2, Ter2); legend TOGGLE

Tr3=feedback(1.3*CrG, 1)
Ter3=feedback(1.3*CerG, 1)
step(Td, Tr3, Ter3); legend TOGGLE

Tr4=feedback(1.24*CrG, 1)
Ter4=feedback(1.32*CerG, 1)
step(Td, Tr4, Ter4); legend TOGGLE

stepinfo(Tfreq2,'RiseTimeLimits',[0,1])
stepinfo(Tr4,'RiseTimeLimits',[0,1])
stepinfo(Ter4,'RiseTimeLimits',[0,1])

CTf = feedback(.8*Cfreq,G)
CTr = feedback(1.24*Cr,G)
CTer = feedback(1.32*Cer,G)

step(CTf, CTr, CTer); legend TOGGLE
