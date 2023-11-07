function [TrTime] = TRB(Flow,city)
load(strcat(city,'/Graphs.mat'));
Adj = adjacency(G_road);
Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);
TrTime = G_road.Edges.Weight'.* (1 + 0.15* (Flow'./Capacity).^4 );
end

