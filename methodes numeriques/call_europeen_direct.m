% Set up les parametres
clear all;
clc;
 % Les paramètres
S=100;
risk=0.05;
strike=110;
sigma=0.25;
T=2;
a=0;
b=2*S;

% période 
N=250;
M=250;

% les pas 
h=(b-a)/(N+1);
k=T/(M+1);

%évolution de temps et de l'espace
x=[a:h:b];
t=[0:k:T];

%Construction des matrices
B=zeros(N,N);

for i=1:N-1
    B(i,i)=1+(risk +((sigma*x(i))^2)/(h^2))*k;
    B(i,i+1)= -(((sigma*x(i))^2)/(2*h^2)+ ((risk*x(i)))/(2*h))*k;
    B(i+1,i)=-(((sigma*x(i+1))^2)/(2*h^2)- ((risk*x(i+1)))/(2*h))*k;
   % B(i,1)= sin(pi*h*i),
end 
B(N,N)=1+(risk +((sigma*x(N))^2)/(h^2))*k;


F=zeros(N,M);
 for i=1:M+1
    ti= -risk*t(i+1);

  F(N,i)=((((sigma*x(N))^2)/(2*h*h))+ (risk*x(N))/(2*h))*k*(x(N)-strike*exp(ti));
 end


%matrice U
UU=zeros(N,M);
V=zeros(N+2,M+2);
for i=1:N
    UU(i,1)= max((x(i)-strike),0);
end

for i=1:M+1
    UU(:,i)=UU(:,i)+F(:,i);
    [L,U,sol,er1,er2]=decomp_LU(B,UU(:,i));
    UU(:,i+1)=sol';
end
for i=2:N+1
  for j=1:M+2
   V(i,j)= UU(i-1,j);
  end
end
V(1,:)=0;
for i=1:M+2
     ti= -risk*t(i);
    V(N+2,i)=b-strike*exp(ti);
end


%V la martice solution de l'équation différentielle

w=T/2 ;
s=1;
d= abs(w-t(1));

for i=1:M+2
    if (abs(w-t(i)))<d
      s=i;
      d= abs(w-t(i));
    end
end

w2=S;
s2=1;
d=abs(S-x(1));
for i=1:N+2
    if (abs(S-x(i)))<d
      s2=i;
      d= abs(S-x(i));
    end
end

call=UU(s2+1,s+1)
figure(1);

[call,put]= blsprice(S, strike, risk,T/2,sigma)

surf(t,x,v)
