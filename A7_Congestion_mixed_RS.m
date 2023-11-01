clc
clear all
city='SF';
Delay = 10;
load(strcat(city,'/Graphs.mat'));
Adj = adjacency(G_road);
Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);

if strcmp('SF', city)
DemandS = DemandS/24; %%daily to hourly only for SF!!!!!!!
end

OriginalDemand= full(DemandS);
DemandS = full(DemandS);
TotDems = sum(DemandS,'all');
for mult = [1.2]  %total demands
for PenRate = [0.5 1] %penetration rate
privateBase.Flows = zeros(N_edges,1);
solBase.Flows = zeros(N_edges,1);
counter = 0;
err = 1;
while true
prev = ([privateBase.Flows; solBase.Flows]);
solBase =LTIFM_NRP_congestion(PenRate*mult*OriginalDemand,city,privateBase.Flows);%privateBase.flows);
privateBase = TAP_casadi((1-PenRate)*mult*OriginalDemand,city,solBase.Flows);
counter = counter + 1;
if counter > 1
err = norm([privateBase.Flows; solBase.Flows] - prev)/norm(prev);
end
if err < 1e-2
    counter
    break;
end
end
solBase.Private = privateBase.Flows;
save(strcat(city,'/Congestion/Results/Mixed/solBase_Delay',num2str(Delay),'Exo',num2str(PenRate),'Dem',num2str(mult),'.mat'),'solBase')
end
end

