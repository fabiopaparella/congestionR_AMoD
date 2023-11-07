function [Opti_RP] = LTIFM_RP_congestion(city,Del)

load(strcat(city,'/Graphs.mat'));
load(strcat(city,'/Congestion/L2/MatL2.mat'));

Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);

x = sdpvar(N_edges,N_nodes,'full');
x_r = sdpvar(N_edges,1,'full');
gamma = sdpvar(N_nodes,N_nodes,N_nodes,N_nodes,3,'full'); % da j a i da j a i
Demm = sdpvar(N_nodes,N_nodes,'full');
D = zeros(N_nodes,N_nodes);
ExogT = sdpvar(N_edges,1,'full'); %0*ones(1,N_edges)';
tau = G_road.Edges.Weight'./ Capacity;
Capacity_th = 1.2*Capacity;
epsilon  = max(0,sum(x,2) + x_r + ExogT - Capacity_th' );


Obj =  G_road.Edges.Weight'*(sum(x,2)+ x_r) + 7.3*tau* ( epsilon.^2 + epsilon.*(Capacity_th' - ExogT ));  %FFT*x
ZeroGamma = zeros(N_nodes,N_nodes,N_nodes,N_nodes,3);
%adding c=0,1,2
for ii1=1:N_nodes
    ii1
    for jj1= 1 : N_nodes
        for ii2=1:N_nodes
            for jj2= 1 : N_nodes
                if ~isempty(sol2_LC{jj1,ii1,jj2,ii2}) & (ii1~=ii2 | jj1~=jj2)
                    if sum(sol2_LC{jj1,ii1,jj2,ii2}(2).Delay < Del) ==2 && sum(sol2_LC{jj1,ii1,jj2,ii2}(3).Delay < Del) ==2
                        D = D + gamma(jj1,ii1,jj2,ii2,1)*sol2_LC{jj1,ii1,jj2,ii2}(1).Dem + ...
                            gamma(jj1,ii1,jj2,ii2,2)*sol2_LC{jj1,ii1,jj2,ii2}(2).Dem + ...
                            gamma(jj1,ii1,jj2,ii2,3)*sol2_LC{jj1,ii1,jj2,ii2}(3).Dem ;
                        %gamma(jj2,ii2,jj1,ii1,1)*sol2_LC{jj2,ii2,jj1,ii1}(2).Dem + ...
                        %gamma(jj2,ii2,jj1,ii1,2)*sol2_LC{jj2,ii2,jj1,ii1}(3).Dem;
                    else
                        D = D + gamma(jj1,ii1,jj2,ii2,1)*sol2_LC{jj1,ii1,jj2,ii2}(1).Dem;
                        ZeroGamma(jj1,ii1,jj2,ii2,2:3)=1;
                    end

                elseif ii1==ii2 & jj1==jj2
                    Mattt=zeros(N_nodes,N_nodes);
                    Mattt(jj1,jj1)=-1;
                    Mattt(ii1,jj1) =1;
                    D = D + gamma(jj1,ii1,jj2,ii2,3)*Mattt; %gamma(jj1,ii1,jj2,ii2,1)*Mattt + gamma(jj1,ii1,jj2,ii2,2)*Mattt set to zero in next zerogamma
                    ZeroGamma(jj1,ii1,jj2,ii2,1:2)=1;
                else
                    ZeroGamma(jj1,ii1,jj2,ii2,:)=1;
               end
            end
        end
    end
end

for ii1=1:N_nodes
    for jj1= 1 : N_nodes
        for ii2=1:N_nodes
            for jj2= 1 : N_nodes
                if isempty(sol2_LC{jj1,ii1,jj2,ii2})
                    ZeroGamma(jj1,ii1,jj2,ii2,:)=1;
                elseif sum(sol2_LC{jj1,ii1,jj2,ii2}(2).Delay > Del) >= 1 && sum(sol2_LC{jj1,ii1,jj2,ii2}(3).Delay > Del) >=1
                    ZeroGamma(jj1,ii1,jj2,ii2,2:3)=1;
                elseif sum(sol2_LC{jj1,ii1,jj2,ii2}(2).Delay > Del) >= 1
                    ZeroGamma(jj1,ii1,jj2,ii2,2)=1;
                elseif sum(sol2_LC{jj1,ii1,jj2,ii2}(3).Delay > Del) >=1
                    ZeroGamma(jj1,ii1,jj2,ii2,3)=1;
                end  
                if ii1==ii2 & jj1==jj2
                    ZeroGamma(jj1,ii1,jj2,ii2,1:2)=1;
                end
                
                    
            end
        end
    end
end

Cons1 = [ Binc*x == D ]; %% fix
Cons2 = [ x >=  0  ];
Cons3 = [ x_r >= 0 ];
       
Cons4 = [ Binc * ( (sum ( x,2) + x_r))  == 0   ];
Cons5 = [];
for ii1 =1: N_nodes
    for jj1 = 1 : N_nodes
        for ii2 =ii1+1: N_nodes
            for jj2 = jj1+1 : N_nodes
                Cons5 = [Cons5
                    gamma(jj1,ii1,jj2,ii2,1) == 0
                    ];
            end
        end
    end
end
Cons6 = [];
Matt = sum(sum(sum(gamma,5),4),3);
Matt2 = squeeze(sum(sum(sum(gamma,5),2),1));
for ii =1: N_nodes
    for jj = 1 : N_nodes  
        if ii~=jj
        Cons6 = [Cons6
             Matt(jj,ii) + Matt2(jj,ii) == Demm(ii,jj)
             
             ];
        end
    end
end
Cons7 = [gamma >= 0];
Cons8 = [gamma( find(ZeroGamma==1) ) == 0] ;
Cons9 = [ sum(x,2) + x_r + ExogT <= 3*Capacity' ];

Cons =  [Cons1
        Cons2
        Cons3
        Cons4
        Cons5
        Cons6
        Cons7
        Cons8
        Cons9
        ];

options = sdpsettings('verbose',1,'solver','Gurobi', 'showprogress',1);
global Opti_RP

Opti_RP = optimizer(Cons,Obj,options,{Demm,ExogT},{Obj,gamma,x,x_r,D,epsilon});
end

