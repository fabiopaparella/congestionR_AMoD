%% same as A2, but stores every combination, not only the best one
clc
clear all
city = 'SF';
load(strcat(city,'/Graphs.mat'));
load(strcat(city,'/solPart_',city,'.mat'));
Adj = adjacency(G_road);
Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);

%% Layer 2
mkdir(strcat(city,'/Congestion/L2'))
sol2_LC{24,24,24,24}=[]; 
for jj1=1:N_nodes %18
    %NYC
    %counter=1;
    jj1
   for ii1=1:N_nodes%15:N_nodes %1
      for ii2=1:N_nodes %12%
         for jj2=1:N_nodes %20  
            if ~any([ii2==jj2,ii1==jj2,ii2==jj1,ii1==jj1]) %,DemandS(ii1,jj1)==0,DemandS(ii2,jj2)==0
               sol2_LC{jj1,ii1,jj2,ii2} = [LTIFM2_SP_congestion(jj1,ii1,jj2,ii2,solPart)]; % matrix with objective, delays, order
               %sol2_LC{counter+1,:} = [LTIFM2_SP_congestion(jj2,ii2,jj1,ii1,solPart)]; % matrix with objective, delays, order
               %counter = counter+2;
            end  
         end
      end
   end
%sol2_LC( sol2_LC(:,1) == 0,: ) = []; %filter out zero obj
%sol2_LC( sol2_LC(:,2) > 20,: ) = []; %filter out above 20 delay
%sol2_LC( sol2_LC(:,3) > 20,: ) = []; %filter out above 20 delay
   
end

save(strcat(city,'/Congestion/L2/MatL2.mat'),'sol2_LC')