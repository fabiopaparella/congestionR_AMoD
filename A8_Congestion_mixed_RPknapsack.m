clc
clear all
city='SF';
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
%mkdir(strcat(city,'/Results'))
%% Layer 2
ppl = 2;

if ppl ==2
    load(strcat(city,'/MatL2.mat'))
    size2 = size(sol2_LC,1);
    Sol2 = [sol2_LC(:,1:3),NaN(size2,2),sol2_LC(:,4:7),NaN(size2,4),sol2_LC(:,8:11),NaN(size2,4)];
    FullList = [Sol2];
end

FullList = sortrows(FullList,1);

for iiii = 1:size(FullList,1)
    vect =  FullList(iiii,6:13)';
    vectR = vect(~isnan( vect));
    vectR = reshape(vectR,2,[])';
    num = size(vectR,1);
    for iii=1:num
    if DemandS(vectR(iii,2),vectR(iii,1)) == 0
    FullList(iiii,1) = 0;
    end
    end
end
FullList(FullList(:,1) >= -0.01,:) = [];
%%
CountGamma = [];
for WaitingTime = [20] %in min 2 5 10 15
    for Delay = [10] % in min
        for PenRate= [0.5:0.1:1]
            for mult= [1 1.5 2] %0.0078 0.0156 0.0312 0.0625 0.125 0.25 0.5 1 2
                gamma = 0;
                TotGamma2 = 0;
                TotGamma3 = 0;
                TotGamma4 = 0;
                Cumul_delay2 = 0;
                Cumul_delay3 = 0;
                Cumul_delay4 = 0;
                DemandS =  PenRate*mult* OriginalDemand;
                Demands_rp = zeros(N_nodes,N_nodes);
                
                for iii = 1:size(FullList,1)
                    if isnan(gamma)
                        iii
                        break
                    end
                    if isnan(FullList(iii,4)) && isnan(FullList(iii,5)) && FullList(iii,2) < Delay && FullList(iii,3) < Delay && DemandS(FullList(iii,7),FullList(iii,6)) >= 10e-5 && DemandS(FullList(iii,9),FullList(iii,8)) >= 10e-5
                        
                        jj1 = FullList(iii,6);ii1 = FullList(iii,7);jj2 = FullList(iii,8); ii2 = FullList(iii,9);
                        gamma = min([DemandS(ii1,jj1),DemandS(ii2,jj2)])*probcombN([DemandS(ii1,jj1),DemandS(ii2,jj2)],WaitingTime)/2;
                        
                        Gamma0 = zeros(N_nodes,N_nodes);
                        if FullList(iii,14:17) == [1 2 1 2]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(ii1,jj2) = 1;
                            Gamma0(ii2,ii1) = 1;
                        elseif FullList(iii,14:17) == [1 2 2 1]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(ii2,jj2) = 1;
                            Gamma0(ii1,ii2) = 1;
                        end
                        
                        matrow = [jj1, ii1; jj2, ii2];
                        multip = size(unique(matrow, 'rows', 'first'),1);
                        Demands_rp = Demands_rp +  multip* gamma* Gamma0;
                        DemandS(ii1,jj1) = DemandS(ii1,jj1) - multip*gamma;
                        DemandS(ii2,jj2) = DemandS(ii2,jj2) - multip*gamma;
                        Cumul_delay2 = Cumul_delay2 + multip*gamma* (FullList(iii,2) + FullList(iii,3)) ;
                        TotGamma2 = TotGamma2 + multip*gamma;
                        %CountGamma = [CountGamma gamma];
                    end
                    
                    
                end
                Demands_rp = Demands_rp - diag(diag(Demands_rp));
                privateBase.Flows = zeros(N_edges,1);
                solKN.Flows = zeros(N_edges,1);
                counter = 0;
                err = 1;
                while true
                    prev = ([privateBase.Flows; solKN.Flows]);
                    solKN = LTIFM_NRP_congestion(Demands_rp+DemandS ,city,privateBase.Flows);
                    privateBase = TAP_casadi((1-PenRate)*mult*OriginalDemand,city,solKN.Flows);
                    counter = counter + 1;
                    if counter > 1
                        err = norm([privateBase.Flows; solKN.Flows] - prev)/norm(prev);
                    end
                    if err < 1e-2
                        counter
                        break;
                    end
                    
                end
                NonRP_count = sum(DemandS,'all'); 
                Tot_count = sum(PenRate*mult*OriginalDemand,'all');
                solKN.Private = privateBase.Flows;
                solKN.RatioRP = 1 - NonRP_count/Tot_count;
                save(strcat(city,'/Congestion/Results/Mixed/solKn_Delay',num2str(Delay),'Exo',num2str(PenRate),'Dem',num2str(mult),'.mat'),'solKN')
            end
            
        end
    end
end
