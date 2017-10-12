set V;
set E := {v in V, w in V: v <> w};

param edge_capacity {E} >= 0;

var edge_flow {E} >= 0;

maximize Flow: edge_flow["T", "S"];

subject to Capacity {(v, w) in E}: edge_flow[v, w] <= edge_capacity[v, w];
subject to Balance {v in V}: sum {(v, w) in E} edge_flow[v, w] = sum {(w, v) in E} edge_flow[w, v];

option solver cplex;
