function C = Controle_Kp_Tp(G,Kp,Mp,Tp,T1,T2,T3)
%Mp é o overshoot
%Kp é o Kp desejado que o controlador de avanço possua
%Tp tempo de pico
%T1 retardo antes da malha fechada
%T2 retardo antes da malha malha direta
%T3 retardo antes da malha malha de retorno
fprintf("Não iremos encontrar um controlador lead-lag,\nmas sim um de avanço para o Kp")

razao = evalfr(G,0);
fprintf("G(0) = %f\nlogo ",razao);
razaoKa_b = Kp/razao;
fprintf("k*a/b = %f\n",razaoKa_b);
fprintf("k*a/b = %f\n",razaoKa_b);
fprintf("Encontrando o valor de C\n")
C = Controle_razaoKa_b_Tp(G,razaoKa_b,Mp,Tp,T1,T2,T3);

end