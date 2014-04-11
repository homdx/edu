n = 4000;
A = zeros(4000);
tic
for jt=1:4000
	for it=1:4000
		A(it,jt) = 1;
	end
end
toc