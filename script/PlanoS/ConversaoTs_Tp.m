function [Tp,Mp] = ConversaoTs_Tp(Ts,Mp,x)

xi = -log(Mp)/sqrt(pi^2+[log(Mp)]^2);
wn  = -log(x/100)/((Ts)*xi);

fprintf("\ntem que retornar [Tp,Mp]\nTp = pi/(wn*sqrt(1-xi^2))\n");
Tp = pi/(wn*sqrt(1-xi^2));

end