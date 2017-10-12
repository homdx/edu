var x_B;
var x_C;

maximize Profit: 25 * x_B + 30 * x_C;

s.t. Time: (1/200) * x_B + (1/140) * x_C <= 40;
subject to B_limit: 0 <= x_B <= 6000;
subject to C_limit: 0 <= x_C <= 4000;

option solver cplex;
