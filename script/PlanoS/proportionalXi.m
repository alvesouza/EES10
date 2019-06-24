function K = proportionalXi(G,xi)
[z,p,k] = zpkdata(G);
raizes = cell2mat(z);
polos = cell2mat(p);

j = sqrt(-1);
nroot = size(raizes,1);
npoles = size(polos,1);
q = [-xi+j/sqrt(1-xi^2) 0];

fprintf("\n||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("\nCalculando o G(%s)\n",mat2str(q));
denominador = [1];
i = 1;
while i<= npoles
    
       
    denominador = conv(denominador,conv([q(1) -polos(i)],[real(q(1))-j*imag(q(1)) -(real(polos(i))-j*imag(polos(i)))]));
    i = i + 1;
end

fprintf("\ndenominador = %s\n\n\n",mat2str(denominador));

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
fprintf("\nCalculando o numerador\n");
numerador = [1];
i = 1;
while i<= nroot
    
    fprintf("raiz = %f + i*%f\n",real(raizes(i)),imag(raizes(i)));
    
    numerador = conv(numerador,conv([q(1) -raizes(i)],[real(q(1))-j*imag(q(1)) -(real(polos(i))-j*imag(polos(i)))]));
    i = i + 1;
end

fprintf("\nnumerador = %s + j*\n\n",mat2str(real(numerador)),mat2str(imag(numerador)));

fprintf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n");
wn = roots(imag(numerador));
n = size(wn,1);
i = 1;
k = k*1.^(1:n);
while i<=n
    k(i) = k(i)*polyval(denominador,wn(i))/polyval(numerador,wn(i));
    i = i+1; 
end
fprintf("wn = %s\n",mat2str(wn));
fprintf("K = %s\n",mat2str(k));
K = k;
end