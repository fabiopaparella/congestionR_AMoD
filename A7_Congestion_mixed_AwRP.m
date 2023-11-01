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
privateBase.Flows = zeros(N_edges,1);

[OptiRP] = LTIFM_RP_congestion(city,Delay);

for mult = [1 1.5 2]  %total demands
    for PenRate = [0.5:0.1:1] %penetration rate
        privateBase.Flows = zeros(N_edges,1);
        counter = 0;
        err = 1;
        while true
            prev = ([privateBase.Flows; sum(sol1{3},2) + sol1{4} ]);
            sol1 = OptiRP(PenRate*mult*DemandS,privateBase.Flows);
            privateBase = TAP_casadi((1-PenRate)*mult*OriginalDemand,city, sum(sol1{3},2) + sol1{4});
            counter = counter + 1;
            if counter > 1
                err = norm([privateBase.Flows; solKN.Flows; solBase.Flows] - prev)/norm(prev);
            end
            if err < 1e-2
                counter
                break;
            end
            
        end
        
        solRP.x = sol1{3};
        solRP.xr = sol1{4};
        solRP.obj = sol1{1};
        solRP.Dem = sol1{5};
        solRP.Gamma = sol1{2};
        solRP.eps = sol1{6};
        solRP.Flows = sum(solRP.x,2) + solRP.xr;
        solRP.Private = privateBase.Flows;
        save(strcat(city,'/Congestion/Results/Mixed/solRP_Delay',num2str(Delay),'Exo',num2str(PenRate),'Dem',num2str(mult),'.mat'),'solRP')
        
        
    end
end
