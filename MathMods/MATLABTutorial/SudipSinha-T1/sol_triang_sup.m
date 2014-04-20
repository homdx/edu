function x=sol_triang_sup(U,b)
[n,m]=size(U);
x=zeros(1,n);
x(n)=b(n)/U(n,n);
i=n-1;
while(i>=1)
somma=0;
for k=i+1:n
somma=somma+U(i,k)*x(k);
end;
x(i)=(b(i)-somma)/U(i,i);
i=i-1;
end;
return;
end