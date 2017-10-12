set Products;

param total_hours  >=  0;
param profit_per_mass {Products}  >=  0;
param mass_per_hour {Products}  >=  0;
param capacity {Products}  >=  0;

var mass {Products}  >=  0;

maximize Profit:  sum {p in Products} profit_per_mass[p] * mass[p];

subject to time:  sum {p in Products} mass[p] / mass_per_hour[p]  <=  total_hours;
subject to production_limit {p in Products}:  0  <=  mass[p] <=  capacity[p];

option solver cplex;
