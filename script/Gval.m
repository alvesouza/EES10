function Gval = Gval(G,T1,T2,T3,q)
[z,p,k] = zpkdata(G);
raizes = cell2mat(z);
polos = cell2mat(p);
x = real(q);
y = imag(q);
j = sqrt(-1);
nroot = size(raizes,1);
npoles = size(polos,1);

Gval = 1;

fase = 1;
i = 1;
fprintf("\n||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("\nCalculando o G(%s)\n",mat2str(q));
while i<= npoles
    
    fprintf("polo = %f + i*%f\n",real(polos(i)),imag(polos(i)));
    fprintf("(q-polos(i))=>(%f + j*(%f))\n",real(q-polos(i)),imag(q-polos(i)));
        
    Gval = Gval/(q-polos(i));
    fase = fase*(q-polos(i));
    i = i + 1;
end

fprintf("\n(q-polos) = %f\n\n\n",fase);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("\nCalculando os (q-raizes)\n");
fase = 1;
i = 1;
while i<= nroot
    
    fprintf("raiz = %f + i*%f\n",real(raizes(i)),imag(raizes(i)));
    fprintf("(q-raizes(i))=>(%f + j*(%f))\n",real(q-raizes(i)),imag(q-raizes(i)));

    Gval = Gval*(q-raizes(i));
    fase = fase*(q-raizes(i));
    i = i + 1;
end

fprintf("\n(q-raizes) = %f\n\n",fase);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("abs(k)*exp(-real(q)*(T2+T3))=>abs(%f)*exp((%f+j*%f)*(%f+%f))= %f\n",k,-real(q),imag(q),T2,T3,abs(k)*exp(-real(q)*(T2+T3)));


Gval = Gval*(abs(k)*exp(-q*(T2+T3)));
fprintf("\nCalculando o G(%s) = %f\n",mat2str(q),Gval);
end