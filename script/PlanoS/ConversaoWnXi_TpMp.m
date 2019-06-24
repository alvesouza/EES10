function [Tp,Mp] = ConversaoWnXi_TpMp(wn,xi)
fprintf("\ntem que retornar [Tp,Mp]\n");
Mp = exp(-pi*xi/sqrt(1-xi^2));
Tp = pi/(wn*sqrt(1-xi^2));

end