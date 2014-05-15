%%	Load the data:
load fisheriris

%%	Use the petal length (PL) and petal width (PW) measurements:
PL = meas(:,3);
PW = meas(:,4);

%%	Plot the data, showing the classification:
h1 = gscatter(PL,PW,species,'krb','ov^',[],'off');
set(h1,'LineWidth',2);
legend('Setosa', 'Versicolor', 'Virginica', 'Location', 'best');

%%	Create a linear classifier:
X = [PL,PW];
cls = ClassificationDiscriminant.fit(X,species);

%%	Plot the classification boundaries:
hold on
K = cls.Coeffs(2,3).Const;
L = cls.Coeffs(2,3).Linear;
% Plot the curve K + [x,y]*L = 0:
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2 = ezplot(f,[.9 7.1 0 2.5]);
set(h2,'Color','r','LineWidth',2);

K = cls.Coeffs(1,2).Const;
L = cls.Coeffs(1,2).Linear;
% Plot the curve K + [x1,x2]*L = 0:
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h3 = ezplot(f,[.9 7.1 0 2.5]);
set(h3,'Color','k','LineWidth',2);
axis([.9 7.1 0 2.5]);
xlabel('Petal Length');  ylabel('Petal Width');
title('{\bf Linear Classification with Fisher Training Data}');
