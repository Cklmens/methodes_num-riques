% Set up les parametres
a=0   % valeur initial
b=1
T=1   % temps
c=1/(pi*pi) 

% période 
N=100
M=1000

% les pas 
h=(b-a)/(N+1)
k=T/(M+1)

%évolution de temps et de l'espace
x=[a,h,b]
t=[0,k,T]

