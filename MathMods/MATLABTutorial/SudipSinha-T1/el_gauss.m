function [U,b]=el_gauss(A,b)
A=[A,b'];
[n,m]=size(A);
for k=1:n-1
for i=k+1:n
Mik=A(i,k)/A(k,k);
for j=k:m
A(i,j)=A(i,j)-Mik*A(k,j);
end;
end;
end;
U=A(:,1:m-1);
c=A(:,m);
b=c';
end