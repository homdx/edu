year = [1975 1980 1985 1990];
west = [72.8 74.2 75.2 76.4];
east = [70.2 70.2 70.3 71.2];

yeari = [1970 1975 1980 1983 1985 1988 1990 1995];
westil = interp1(year, west, yeari, 'linear', 'extrap');
eastil = interp1(year, east, yeari, 'linear', 'extrap');
westis = spline(year, west, yeari);
eastis = spline(year, east, yeari);

subplot(1, 2, 1);
plot(yeari, westil, yeari, eastil);
title('Linear'); xlabel('Year'); ylabel('Life expectancy (years)');
legend('West', 'East'); grid on;

subplot(1, 2, 2);
plot(yeari, westis, yeari, eastis);
title('Spline'); xlabel('Year'); ylabel('Life expectancy (years)');
legend('West', 'East'); grid on;
