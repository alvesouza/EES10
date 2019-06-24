function Gval = Gval(G,T1,T2,T3,s)
[z,p,k] = zpkdata(G);
raizes = cell2mat(z);
polos = cell2mat(p);
x = real(s);
y = imag(s);
j = sqrt(-1);
nroot = size(raizes,1);
npoles = size(polos,1);

Gval = 1;

fase = 1;
i = 1;
fprintf("\n||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("\nCalculando o G(%s)\n",mat2str(s));
while i<= npoles
    
    fprintf("polo = %f + i*%f\n",real(polos(i)),imag(polos(i)));
    fprintf("(s-polos(i))=>(%f + j*(%f))\n",real(s-polos(i)),imag(s-polos(i)));
        
    Gval = Gval/(s-polos(i));
    fase = fase*(s-polos(i));
    i = i + 1;
end

fprintf("\n(s-polos) = %f\n\n\n",fase);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("\nCalculando os (s-raizes)\n");
fase = 1;
i = 1;
while i<= nroot
    
    fprintf("raiz = %f + i*%f\n",real(raizes(i)),imag(raizes(i)));
    fprintf("(q-raizes(i))=>(%f + j*(%f))\n",real(s-raizes(i)),imag(s-raizes(i)));

    Gval = Gval*(s-raizes(i));
    fase = fase*(s-raizes(i));
    i = i + 1;
end

fprintf("\n(q-raizes) = %f\n\n",fase);

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("(k)*exp(-real(s)*(T2+T3))=>abs(%f)*exp((%f+j*%f)*(%f+%f))= %f\n",k,-real(s),imag(s),T2,T3,(k)*exp(-real(s)*(T2+T3)));


Gval = Gval*((k)*exp(-s*(T2+T3)));
fprintf("\nCalculando o G(%s) = %f\n",mat2str(s),Gval);
end