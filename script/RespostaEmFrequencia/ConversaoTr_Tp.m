function [Tp,Mp] = ConversaoTr_Tp(Tr,Mp)

xi = -log(Mp)/sqrt(pi^2+[log(Mp)]^2);
wn  = (pi-acos(xi))/((Tr)*sqrt(1-xi^2));

fprintf("\ntem que retornar [Tp,Mp]\nfc = 2*pi/");
Mp = exp(-pi*xi/sqrt(1-xi^2));
Tp = pi/(wn*sqrt(1-xi^2));

end