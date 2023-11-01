clc
clear all
city='NYC120'; %% NYC120 or SF
load(strcat(city,'/Graphs.mat'));
Adj = adjacency(G_road);
Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);
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

%% casadi optimizer

import casadi.*
clear opti
opti = casadi.Opti()

x = opti.variable(N_edges,N_nodes);
y = opti.variable(N_edges,1);
ExogT_par = opti.parameter(N_edges,1);
Demands_par = opti.parameter(N_nodes,N_nodes);

Obj = G_road.Edges.Weight'*( y + Capacity'.^(-4)*0.03.*( (y + ExogT_par).^5 - ExogT_par.^5)) ;
Cons = [];
Cons1 = [{ Binc*x == Demands_par }];
Cons2 = [{ x(:) >=  0 }];
Cons3 = [{ sum(x,2) + ExogT_par <= 3*Capacity' }];
Cons4 = [{sum(x,2) == y}];
Cons = [Cons1;
        Cons2;
        Cons3;
        Cons4;
          ];
opti.minimize(Obj);
opti.subject_to(Cons);
solver_options.expand = 1;
solver_options.print_time = 1;
solver_options.ipopt.linear_solver = 'mumps';
opti.solver('ipopt', solver_options);
%%
CountGamma = [];
for WaitingTime = [15] %in min 2 5 10 15
    for Delay = [10] % in min
        for PenRate= [0.2] % 0.01 0.1:0.1:0.9 0.99
            for Psi = [0.2]%0.01 0.1:0.1:0.9 0.99
            for mult= [4] %0.0078 0.0156 0.0312 0.0625 0.125 0.25 0.5 1 2
                gamma = 0;
                TotGamma2 = 0;
                TotGamma3 = 0;
                TotGamma4 = 0;
                Cumul_delay2 = 0;
                Cumul_delay3 = 0;
                Cumul_delay4 = 0;
                DemandS =  PenRate*Psi*mult*OriginalDemand;
                Demands_rp = zeros(N_nodes,N_nodes);
                Delay_cong =  FullList(:,1);
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
                        Delay_cong(iii,1) = gamma;
                    end
                    
                    
                end
                
                Demands_rp = Demands_rp - diag(diag(Demands_rp));
                privateBase.Flows = zeros(N_edges,1);
                solKN.Flows = zeros(N_edges,1);
                solBase.Flows = zeros(N_edges,1);
                counter = 0;
                err = 1;
                priv_demands = OriginalDemand;
                for ii=1:N_nodes
                priv_demands(ii,ii) = -sum(priv_demands(:,ii)) - priv_demands(ii,ii);
                end
                
                while true
                    prev = ([privateBase.Flows; solKN.Flows ; solBase.Flows]);
                    opti.set_value(Demands_par, (1-PenRate)*mult*priv_demands);
                    opti.set_value(ExogT_par, solKN.Flows+solBase.Flows);
                    privateBase1 = opti.solve();
                    privateBase.x = privateBase1.value(x);
                    privateBase.Flows = privateBase1.value(y);
                    privateBase.obj = privateBase1.value(Obj);
                    solKN = LTIFM_NRP_congestion(Demands_rp+DemandS,city,privateBase.Flows+solBase.Flows);
                    solBase = LTIFM_NRP_congestion(PenRate*(1-Psi)*mult*OriginalDemand,city,privateBase.Flows+solKN.Flows);
                    counter = counter + 1;
                    if counter > 1
                        err = norm([privateBase.Flows; solKN.Flows; solBase.Flows] - prev)/norm(prev);
                        counter
                    end
                    if err < 1e-2 | counter > 10
                        counter
                        break;
                    end
                    
                end
                
                NonRP_count = sum(DemandS,'all'); 
                Tot_count = sum(PenRate*Psi*mult*OriginalDemand,'all');
                solKN.RatioRP = 1 - NonRP_count/Tot_count;
                solKN.gammas = Delay_cong;
                save(strcat(city,'/Congestion/Results/3Mixed/Delay',num2str(Delay),'Exo',num2str(PenRate),'RP',num2str(Psi),'Dem',num2str(mult),'.mat'),'solKN','privateBase' ,'solBase')
            end
            end  
        end
    end
end
