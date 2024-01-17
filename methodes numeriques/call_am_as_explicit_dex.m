% Set up les parametres
clear all;
clc;
 % Les paramètres
S=100;
Smin=S-S/2;
Smax=S+S/2;
risk=0.02;
strike=100;
sigma=0.3;
T=2;

%a=0;
%b=2*S;
a=-100;
b=-0;

% période 
N=500;
M=2.1*N*N;

% les pas 
h=(b-a)/(N+1);
k=T/(M+1);
ho=(-Smin+Smax)/(N+1);
So=[Smin:ho:Smax];

% Inita=ialiser x
%x=zeros(N+1);
%évolution de temps et de l'espace
for i=1:N+2
   % x(i)=-(mean(So(1:i)))/So(i);
end

x=[a:h:b];
t=[0:k:T];

%Construction des matrices
B=zeros(N,N);


F=zeros(N,M);
 for i=1:M+1
  F(N,i)=((((sigma*x(N))^2)/(2*h*h))- (risk*x(N))/(2*h) - (1+x(N))/(2*h*(T-t(i))))*k*x(N);
 end


%matrice U
UU=zeros(N,M);
V=zeros(N+2,M+2);
O=0;
for i=1:N
    UU(i,1)= max((-x(i)-1),0);
end

for j=1:M+1

 for i=1:N-1
    B(i,i)=1-(((sigma*x(i))^2)/(h^2)- abs(((risk*x(i)))/(2*h)+ (1+x(i))/(2*h*(T-t(j)))))*k;
    B(i,i+1)= (((sigma*x(i))^2)/(2*h^2))*k;
    B(i+1,i)=(((sigma*x(i+1))^2)/(2*h^2)- abs(((risk*x(i+1)))/(2*h)+(1+x(i+1))/(2*h*(T-t(j)))))*k;
   % B(i,1)= sin(pi*h*i),
end 
B(N,N)=1-(((sigma*x(N))^2)/(h^2)+ ((risk*x(N)))/(2*h)+ (1+x(N))/(2*h*(T-t(N))))*k;

 UU(:,j+1)= B*UU(:,j)+F(:,j) ;
O=O+1;
   % si dans le cas américain on ajoute:
    for i=1:N
        man= max(-x(i)-1,0);
        UU(i,j+1)=max(UU(i,j+1),man);
    end

end

s2=1;
d= abs(-1-x(1));

for i=1:N+1
    if (abs(-1-x(i)))<d
      s2=i;
      d= abs(-1-x(i));
    end
end
w=T/2 ;
s=1;
d= abs(w-t(1));

for i=1:M+2
    if (abs(w-t(i)))<d
      s=i;
      d= abs(w-t(i));
    end
end

CALL=S*UU(s2,s)