%o rlocus ele serve para uma dada eq caract.
%assim, precisamos encontrar algum G tal que possua
    % o mesmo rlocus, ou seja, o mesmo den 
    % em malha fechada
%E que possua o parametro que desejamos variar   
%Exemplo:

%variando b:

b = 0.263
G9 = zpk([-10],[0 -b -8],[10])

G9_ = tf([1 8 0], [1 8 10 100])

rlocus(G9_)

roots([-j -(8+b) j*(10+8*b) 100])

%varianda a:

a = 1
G9_ = tf([10],[1 9 18 0])

rlocus(G9_)

