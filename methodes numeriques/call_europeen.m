% Set up les parametres
clear all;
clc;
 % Les paramètres
S=100;
risk=0.05;
strike=110;
sigma=0.25;
Smin=S-S/2;
Smax=S+S/2;
a= log(Smin/strike);  % valeur initial
b=log(Smax/strike);
q=(2*risk)/(sigma^2);
T=2;   % temps
lo=0;
lf=(sigma^2)*(T/2);
t1=T/2;
c=1 ;

% période 
N=250;
M=250;

% les pas 
h=(b-a)/(N+1);
k=(lf-lo)/(M+1);

%évolution de temps et de l'espace
x=[a:h:b];
l=[lo:k:lf];

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
F=zeros(N,M);
V=zeros(N+2,M+2);
for i=1:N
    temp=0.5*(q-1)*x(i);
    UU(i,1)=exp(temp)*max((exp(x(i))-1),0);
end

for i=1:M+1
   tmp1=-risk*2*l(i)/(sigma^2);
tmp2=0.5*(q-1)*b+0.25*((q+1)^2)*l(i);
F(N,i)=(exp(b)-exp(tmp1))*exp(tmp2)*c*k*(1/h^2); % ajouter 1
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
   tmp1=-2*risk*l(i)/(sigma^2);
tmp2=0.5*(q-1)*b+0.25*((q+1)^2)*l(i);
V(N+2,i)=(exp(b)-exp(tmp1))*exp(tmp2);
end

%V la martice solution de l'équation différentielle

w=((sigma^2)*(T-t1))/2 ;
s=1;
d= abs(w-l(1));

for i=1:M+2
    if (abs(w-l(i)))<d
      s=i;
      d= abs(w-l(i));
    end
end

w2=log(S/strike);
s2=1;
d=abs(w2-x(1));
for i=1:N+2
    if (abs(w2-x(i)))<d
      s2=i;
      d= abs(w2-x(i));
    end
end



W=zeros(N+2,M+2);

for i=1:N+2
    for j=1:M+2
        tmp= 0.5*(q-1)*x(i)+0.25*((q+1)^2)*l(j);
        W(i,j)=V(i,j)*strike*exp(-tmp);

end
end

tmpe= 0.5*(q-1)*w2+0.25*((q+1)^2)*w;
call=V(s2+1,s+1)*strike*exp(-tmpe)


figure(1);

[call,put]= blsprice(S, strike, risk,T/2 ,sigma)
W(s2-2,s-2);
surf(l,x,W)
