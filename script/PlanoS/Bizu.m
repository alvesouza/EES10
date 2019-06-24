%%
%nyquist
%zeros 1+GH(Z)[polos da equação de interesse]
%polos 1+GH(p)
%numero de circulações anti-horarias em torno de -1/k
%P = N+Z
%%
%ressonancia
%Asaida = Ainp*|G(jw)|
%phi = phase(G(jw))(atraso)
%%
%ressonancia
%Mp = 1/(2*xisqrt(1-xi^2))%margem de pico(não é em db)
%wp = wn*sqrt(1-2*xi^2)%frequencia de pico(não é em rad/s)
%wbw = wn*sqrt((1-2*xi^2)+sqrt(4*xi^4-4xi^2+2))
%%
%PI
%Improve steady-state error
%1. Increases system type.
%2. Error becomes zero.
%3. Zero at -zc is small and negative.
%4. Active circuits are required to implement.

%Lag
%Improve steady-state error
%1. Error is improved but not driven to zero.
%2. Pole at -pc is small and negative.
%3. Zero at -zc is close to, and to the left of, the
%pole at -pc.
%4. Active circuits are not required to implement.

%PD
%Improve transient response
%1. Zero at zc is selected to put design point on
%root locus.
%2. Active circuits are required to implement.
%3. Can cause noise and saturation; implement
%with rate feedback or with a pole (lead).

%Lead
%Improve transient response
%1. Zero at -zc and pole at -pc are selected to put
%design point on root locus.
%2. Pole at -pc is more negative than zero at -zc.
%3. Active circuits are not required to implement.

%PID
%Improve steady-state error and transient response
%1. Lag zero at -zlag and pole at origin improve steady-state error.
%2. Lead zero at -zlead improves transient response.
%3. Lag zero at -zlag is close to, and to the left of, the origin.
%4. Lead zero at -zlead is selected to put design point on root locus.
%5. Active circuits required to implement.
%6. Can cause noise and saturation; implement with rate feedback or with an additional pole.

%Lag-lead
%Improve steady-state error and transient response
%1. Lag pole at -plag and lag zero at -zlag are used to improve steady-state error.
%2. Lead pole at -plead and lead zero at -zlead are used to improve transient response.
%3. Lag pole at -plag is small and negative.
%4. Lag zero at -zlag is close to, and to the left of, lag pole at -plag.
%5. Lead zero at -zlead and lead pole at -plead are selected to put design point on root locus.
%6. Lead pole at -plead is more negative than lead zero at -zlead.
%7. Active circuits are not required to implement.
%%
nichols(sys)
nichols(sys,w)
controlSystemDesigner(G)
mat2str([0 0 0;0 0 0])
polyval(p,2)%%p é vetor de um polinomio e x = 2
conv(p1,p2)%%multiplica polinomios
roots(p)
inv(A)%%inverte a matriz
sys = tf(num,dem)
bode(sys)
nyquist(sys)
T = feedback(G*k,1);
step(sys)
step(sys,Tfinal)
evalfr(G*C,q)
%%
%resolver equações
%A*X = B
%X = inv(A)*B

%%formulas
r = -xi + i*sqrt(1-xi^2)
xi = sqrt(log(Mp)^2/(log(Mp)*log(Mp) + pi*pi))
pol = -wn*xi + wn*sqrt(1-xi^2)*i
Mp = exp(-xi*pi/sqrt(1-xi^2))
stepinfo(T,'RiseTimeLimits',[0,1])
[z,p,k] = zpkdata(C*G2);
z = cell2mat(z);
p = cell2mat(p);