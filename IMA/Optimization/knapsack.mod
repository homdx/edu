set Items;

param capacity >= 0;
param mass {Items} >= 0;
param benefit {Items} >= 0;

var take {Items} binary;

maximize Benefit: sum {i in Items} benefit[i] * take[i];
subject to Capacity: sum {i in Items} mass[i] * take[i]  <=  capacity;

option solver cplex;
