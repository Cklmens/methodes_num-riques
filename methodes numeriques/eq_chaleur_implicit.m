% Set up les parametres
clear all;
clc;

a=0 ;  % valeur initial
b=1;
T=1;   % temps
c=1/(pi*pi) ;

% période 
N=10;
M=10;

% les pas 
h=(b-a)/(N+1);
k=T/(M+1);

%évolution de temps et de l'espace
x=[a:h:b];
t=[0:k:T];

%Construction des matrices
B=zeros(N,N);

for i=1:N-1
    B(i,i)=1+2*(k*c)/(h*h);
    B(i,i+1)= -(k*c)/(h*h);
    B(i+1,i)= -(k*c)/(h*h);
   % B(i,1)= sin(pi*h*i),
end 
B(N,N)=1+2*(k*c)/(h*h);

%matrice U
UU=zeros(N,M);
V=zeros(N+2,M+2);
for i=1:N
    UU(i,1)= sin(pi*(a+h*i));
end

for i=1:M+1
    [L,U,sol,er1,er2]=decomp_LU(B,UU(:,i));
    UU(:,i+1)=sol';
end

for i=2:N+1
  for j=1:M+2
   V(i,j)= UU(i-1,j);
  end
end
V(1,:)=0;
V(N+2,:)=0;
V
for i=1:N+2
    for j=1:M+2
        uexact(i,j)=exp(-((j-1)*k))*sin((a+(i-1)*h)*pi);
    end
end

figure(1);
surf(t,x,V)
figure(2);
surf(t,x,uexact)