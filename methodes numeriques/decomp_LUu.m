function [L,U,sol,err1,err2] = decomp_LUu(A)
%{ n=size(A);
U=zeros(n);
L=eye(n);
U(1,:)=A(1,:);
for i=2:n(1)
    L(i,1)=A(i,1)/U(1,1);
end

for i=2:n(1)
    for j=i:n(2)
     L(i,j-1)=(A(i,j-1)) - L(i,1:i-1)*U(1:i-1,j-1)/U(j-1,j-1)
     U(i,j)=A(i,j) -L(i,:)*U(:,j)
    
    end
end
%}