close all
clear all
clc


Mp=.11; 
Tp=.67;

xi=-log(Mp)/sqrt(pi^2+[log(Mp)]^2)
wn  = (pi)/Tp/sqrt(1-xi^2)

PMd = 100*xi
wc = wn

q = -xi*wn + j*wn*sqrt(1-xi^2)

G = zpk([],[-1 -2],[5])
F = zpk([],[],[1])
H = zpk([],[],[1])

modGq = 5/(abs(q+1)*abs(q+2))
phaseGq = -phase(q+1) - phase(q+2)

s = tf('s')
tau = 0.2
delay = exp(-tau*s)
delay_pade = pade(delay,4)
H = series(H,delay_pade)

Aux1 = series(G,F)

rlocus(Aux1)

zero = -1
pole = -2 + ((-xi*wn - (-2))*2)

k = abs(q-pole)/(modGq*abs(q-zero))

C = zpk([zero],[pole],[k])

%Aux1 = series(C,Aux1)

%T = feedback(Aux1,H)

%stepinfo(T)

angulop = 180 + phaseGq*180/pi + phase(q+zero)*180/pi - tau*imag(q)*180/pi

pole_delay = real(q) - (imag(q)/tand(angulop))

k_delay = exp(tau*real(q))/(modGq)

C_delay = zpk([zero],[pole_delay],[k_delay])

Aux1 = series(C_delay,Aux1)

rlocus(Aux1*delay_pade)

T = feedback(Aux1,H)

stepinfo(T)


