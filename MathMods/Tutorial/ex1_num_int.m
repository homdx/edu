clc;
clear all;
funz=inline('exp(x)');
a=0;
b=10;
tol=0.1;
htr=sqrt((12*tol)/(10*exp(10)));
hsi=2*(((180*tol)/(10*exp(10)))^(0.25));
ntr=round((b-a)/htr);
nsi=round((b-a)/hsi);
fprintf('\n');
fprintf('______________________________________________________________\n');
fprintf('N   \t I_tr      \t\t  Err_tr \t\t\t|\n');
fprintf('______________________________________________________________\t|\n');
I=22025.466; 
tic
for n=1:ntr
 [Itr]=ctrapezi(funz,a,b,n);
 Etr=Itr-I;
fprintf('%i \t| %14.6f \t | %14.6f \t\t|\n',n,Itr,Etr);
end
toc
fprintf('\n');
fprintf('______________________________________________________\n');
fprintf('N   \t I_si      \t\t  Err_si \t\t|\n');
fprintf('______________________________________________________\t|\n');
I=22025.466;  
tic
for n=1:nsi
  [Isi] = simpsonc(funz,a,b,n);
  Esi=Isi-I;
fprintf('%i \t| %14.6f \t | %14.6f \t|\n',n,Isi,Esi);
end
toc