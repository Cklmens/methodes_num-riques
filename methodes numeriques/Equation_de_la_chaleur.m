% Set up les parametres
a=0 ;  % valeur initial
b=1;
T=1;   % temps
c=1/(pi*pi) ;

% période 
N=10;
M=100;

% les pas 
h=(b-a)/(N+1);
k=T/(M+1);

%évolution de temps et de l'espace
x=[a:h:b];
t=[0:k:T];

%Construction des matrices
B=zeros(N+1,N+1);

for i=1:N
    B(i,i)=1-2*(k*c)/(h*h);
    B(i,i+1)= (k*c)/(h*h);
    B(i+1,i)= (k*c)/(h*h);
   % B(i,1)= sin(pi*h*i),
end 
B(N+1,N+1)=1-2*(k*c)/(h*h);

%matrice U
U=zeros(N+1,M+1);
V=zeros(N+2,M+2);
for i=1:N+1
    U(i,1)= sin(pi*h*i);
end

for j=2:M+2
   U(:,j)= B*U(:,j-1);

end
%recopier dans V

for i=1:N+1
  for j=1:M+2
   V(i+1,j)= U(i,j);
  end
end

figure(1)
surf(t,x,V)



