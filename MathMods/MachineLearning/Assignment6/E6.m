X=[-3 3;-3 2;-2 3;-1 1;1 3;2 2;2 3;3 1];
Y=[-1 -1 -1 -1 1 1 1 1]';

C=1;

[n,d]=size(X);

cvx_begin quiet
    variables w(d) b ksi(n)
    dual variable lambda
    minimize 1/2 * sum(w .* w) + C/n * sum(ksi)
    lambda : Y .* (X * w + b) >= 1 - ksi;
    ksi >= 0;
cvx_end

disp('SVM (indices):');
find(lambda' > 1e-6)

if ((sum((lambda .* (Y .* (X*w+b) - 1 + ksi))' < 1e-6)) == length(lambda))
    disp('KKT condition is satisfied');
end

K=X*X';
cvx_begin quiet
    variables alph(n)
    maximize sum(alph) - 0.5 * quad_form(Y .* alph, K)
        alph >= 0;
        alph <= C/n;
        alph'*Y == 0;
cvx_end

if (sum(abs(lambda'-alph')<1e-6)==length(lambda))
    disp('Alpha and Lambda are almost equal');
end

w2 = [sum(alph .* Y .* X(:,1)); sum(alph .* Y .* X(:,2))];

if (abs(w-w2) <= 1e-5)
    disp('Values for w found from cvx and from the formula are almost equal')
end

ac = find(lambda > 1e-6);
b2 = 1/Y(ac(1)) - sum(Y .* alph .* (X * (X(ac(1), :))'));

if (abs(b-b2) <= 1e-4)
    disp('Values for b found from cvx and from the formula are almost equal')
end
