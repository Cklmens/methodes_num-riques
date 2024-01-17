function [L,U,sol,er1,er2]=decomp_LU(A,b)


%le but est de décomposer une matrice A sous la forme L*U
[n,m]=size(A); A0=A; L=eye(n);

for k=1:n-1
  
   M=eye(n);
   N=eye(n);
   for i=k+1:n
      M(i,k)=-A(i,k)/A(k,k);
      N(i,k)= -M(i,k);
   end
   A=M*A;
  
   L=L*N;
end
U=A;

er1=norm(A0-L*U,'fro');
%[y]=descente(L,b);
y(1)=b(1);
for i=2:n
    s=0;
    for j=1:i-1
        s=s+L(i,j)*y(j);
    end
    y(i)=b(i)-s;
    
end

%[sol]=remonte(U,y);
sol(n)=y(n)/U(n,n);
for i=n-1:-1:1
    s=0;
    for j=i+1:n
        s=s+U(i,j)*sol(j);
    end
    sol(i)=(y(i)-s)/U(i,i);
    
end
%er2=norm(b-A0*sol');
L;
U;
sol;
er1;
er2=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

