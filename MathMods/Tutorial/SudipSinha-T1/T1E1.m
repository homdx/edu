%	Test 1 - Exercise 1
clc;
clear;

A = [1 4 3 2 9 5 1 6 5 4;
	 0 1 0 3 0 0 5 0 6 0;
	 2 3 1 0 3 0 0 5 0 6;
	 0 2 0 1 0 3 0 0 5 0;
	 0 9 2 1 1 0 3 0 0 5;
	 0 6 0 2 0 1 0 3 0 0;
	 4 1 0 1 2 0 1 0 3 0;
	 0 4 0 1 0 2 0 1 0 3;
	 1 1 4 1 0 0 2 0 1 0;
	 1 1 4 4 0 0 0 2 0 1];
 
%	Part 1
b = A(1,:);
c = A(:,2);

%	Part 2
sb = sumex1(b);
sc = sumex1(c);
pb = prodex1(b);
pc = prodex1(c);

%	Part 3
b
c

fprintf('Sum of b = %i\n', sb);
fprintf('Sum of c = %i\n', sc);

fprintf('Product of b = %i\n', pb);
fprintf('Product of c = %i\n', pc);
