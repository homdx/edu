set Locations;
set Stores;

param opening_cost {Locations} >= 0;
param demand {Stores} >= 0;
param transport_cost {Locations, Stores} >= 0;

var open {Locations} binary;
var amount_location_to_store {Locations, Stores} >= 0;    # Amount of products shipped from Location l to Store s.

minimize Costs:  
	sum {l in Locations} opening_cost[l] * open[l] +
	sum {l in Locations, s in Stores} transport_cost[l, s] * amount_location_to_store[l, s];

subject to Demands {s in Stores}: sum {l in Locations} open[l] * amount_location_to_store[l, s]  >=  demand[s];
subject to Bounds {l in Locations, s in Stores}: amount_location_to_store[l, s]  <=  demand[s] * open[l];

option solver cplex;
