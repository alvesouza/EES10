function C = PontosInteresseNyquist(G)

j = sqrt(-1);

%separando a função G
[z,p,k] = zpkdata(G);
raizes = cell2mat(z);
polos = cell2mat(p);
fprintf("Determinando os valores de G\n");
fprintf("k = %d\n",k);
fprintf("raizes = %s\n",mat2str(raizes));
fprintf("polos = %s\n\n",mat2str(polos));
nroot = size(raizes,1);
npoles = size(polos,1);
%encontrando o numerador e a derivada do numerador
fprintf("Multiplicando por dem' em cima e embaixo.\n");
numerador = [1];
i = 1;
while i<= nroot
    numerador = conv(numerador,[j -raizes(i)]);
    i = i + 1;
end

%multiplicando o dem' com o num
i = 1;
while i<= npoles
    numerador = conv(numerador,[-j -polos(i)]);
    i = i + 1;
end
denominador = [1];
i = 1;
while i<= npoles
    denominador = conv(denominador,[1 0 (polos(i))^2]);
    i = i + 1;
end
fprintf("Temos G(j*w):\n");
i = 1;
n = size(denominador,1);
fprintf("(%f/(",k);
while i<=n
    if i>1
        fprintf("+");
    end
    if n-i >0
        fprintf("(%f+j*(%f))*s^%d",real(denominador(i)),imag(denominador(i)),(n-i))
    else
        fprintf("(%f+j*(%f))",real(denominador(i)),imag(denominador(i)))
    end
    i = i+1;
end

fprintf("))*A\n\n")
fprintf("//////////////////////////////////////////////\n\n")
fprintf("Parte do numerador com imaginario num*dem':\n")
fprintf("A = ")
i = 1;
n = size(numerador,1);
while i<=n
    if i>1 
        if real(numerador(i))> 0
            fprintf(" + ");
        elseif real(numerador(i))< 0
            fprintf(" - ");
        end
    end
    if n-i >0
        fprintf("%f*s^%d",abs(real(numerador(i))),(n-i))
    else
        fprintf("%f",abs(real(numerador(i))))
    end
    i = i+1;
end
fprintf(" + j*(");
%Parte imaginaria
i = 1;
n = size(numerador,1);
while i<=n
    if i>1 
        if imag(numerador(i))> 0
            fprintf(" + ");
        elseif imag(numerador(i))< 0
            fprintf(" - ");
        end
    end
    if n-i >0
        fprintf("%f*s^%d",abs(imag(numerador(i))),(n-i))
    else
        fprintf("%f",abs(imag(numerador(i))))
    end
    i = i+1;
end
fprintf(")\n");

fprintf("//////////////////////////////////////////////\n\n")
fprintf("Encontrando os w quando o G(jw) for real\n")
fprintf("Im(A) = 0\n")

emReal = roots(imag(numerador));
fprintf("Seus pontos:\n");
i = 1;
n = size(emReal,1);
while i<=n
    A = k*polyval(numerador,emReal(i))/polyval(denominador,emReal(i));
    fprintf("(w=%f,G(jw)=%f+j*(%f));",emReal(i),real(A),imag(A));
    i = i+1;
end

fprintf("//////////////////////////////////////////////\n\n")
fprintf("Encontrando os w quando o G(jw) for imaginario\n")
fprintf("Real(A) = 0\n")

emImag = roots(imag(numerador));
fprintf("Seus pontos:\n");
i = 1;
n = size(emImag,1);
while i<=n
    B = k*polyval(numerador,emImag(i))/polyval(denominador,emImag(i));
    fprintf("(w=%f,G(jw)=%f+j*(%f));",emImag(i),real(B),imag(B));
    i = i+1;
end
fprintf("/////////////////////////////////////////\n")
fprintf("Desenhando o nyquist em ordem crescente com w positivo,com zero e infinito\n")
fprintf("Caso passe por algum polo,ele tem que contornar infinitesimalmente\n")
fprintf("\n")
end