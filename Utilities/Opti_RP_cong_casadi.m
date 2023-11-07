import casadi.*
opti = casadi.Opti()
city = 'SF'
load(strcat(city,'/Graphs.mat'));
load(strcat(city,'/Congestion/L2/MatL2.mat'));
Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);

x = opti.variable(N_edges,N_nodes);
x_r = opti.variable(N_edges,1);
gamma1 = opti.variable(N_nodes,N_nodes,N_nodes,N_nodes); % da j a i da j a i
gamma2 = opti.variable(N_nodes,N_nodes,N_nodes,N_nodes);
gamma3= opti.variable(N_nodes,N_nodes,N_nodes,N_nodes);
Demm = opti.parameter(N_nodes,N_nodes);
D = opti.variable(N_nodes,N_nodes);
Exog_traffic = opti.parameter(N_edges,1); %0*ones(1,N_edges)';
if strcmp('SF', city)
    Capacity = Capacity/24;
else
    Capacity = 2000*ones(1,N_edges);
end
tau = G_road.Edges.Weight'./ Capacity;
Capacity_th = 1.2*Capacity;

epsilon  = max(0,sum(x,2) + x_r + Exog_traffic - Capacity_th' );

Obj =  G_road.Edges.Weight'*(sum(x,2)+ x_r) + tau* ( epsilon.^2 + epsilon.*(Capacity_th' - Exog_traffic ));  %FFT*x
ZeroGamma = zeros(N_nodes,N_nodes,N_nodes,N_nodes,3);
%adding c=0,1,2
for ii1=1:3%N_nodes
    ii1
    for jj1= 1 : N_nodes
        for ii2=1:N_nodes
            for jj2= 1 : N_nodes
                if ~isempty(sol2_LC{jj1,ii1,jj2,ii2}) & (ii1~=ii2 | jj1~=jj2)
                    if sum(sol2_LC{jj1,ii1,jj2,ii2}(2).Delay < Del) ==2 && sum(sol2_LC{jj1,ii1,jj2,ii2}(3).Delay < Del) ==2
                        D = D + gamma1(jj1,ii1,jj2,ii2)*sol2_LC{jj1,ii1,jj2,ii2}(1).Dem + ...
                            gamma2(jj1,ii1,jj2,ii2)*sol2_LC{jj1,ii1,jj2,ii2}(2).Dem + ...
                            gamma3(jj1,ii1,jj2,ii2)*sol2_LC{jj1,ii1,jj2,ii2}(3).Dem ;
                        %gamma(jj2,ii2,jj1,ii1,1)*sol2_LC{jj2,ii2,jj1,ii1}(2).Dem + ...
                        %gamma(jj2,ii2,jj1,ii1,2)*sol2_LC{jj2,ii2,jj1,ii1}(3).Dem;
                    else
                        D = D + gamma1(jj1,ii1,jj2,ii2)*sol2_LC{jj1,ii1,jj2,ii2}(1).Dem;
                        ZeroGamma(jj1,ii1,jj2,ii2,2:3)=1;
                    end

                elseif ii1==ii2 & jj1==jj2
                    Mattt=zeros(N_nodes,N_nodes);
                    Mattt(jj1,jj1)=-1;
                    Mattt(ii1,jj1) =1;
                    D = D + gamma3(jj1,ii1,jj2,ii2)*Mattt; %gamma(jj1,ii1,jj2,ii2,1)*Mattt + gamma(jj1,ii1,jj2,ii2,2)*Mattt set to zero in next zerogamma
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

Cons1 = [ {Binc*x == D} ]; %% fix
Cons2 = [ {x >=  0 } ];
Cons3 = [ {x_r >= 0 }];
       
Cons4 = [ {Binc * ( (sum ( x,2) + x_r))  == 0 }  ];
Cons5 = [];
for ii1 =1: N_nodes
    for jj1 = 1 : N_nodes
        for ii2 =ii1+1: N_nodes
            for jj2 = jj1+1 : N_nodes
                Cons5 = [ {Cons5
                    gamma1(jj1,ii1,jj2,ii2) == 0 }
                    ];
            end
        end
    end
end
Cons6 = [];
Matt = sum(sum(sum(gamma,5),4),3);
Matt2 = squeeze( sum( sum(gamma1+gamma2+gamma3),2),1);
for ii =1: N_nodes
    for jj = 1 : N_nodes  
        if ii~=jj
        Cons6 = [ {Cons6
             Matt(jj,ii) + Matt2(jj,ii) == Demm(ii,jj)
             }
             ];
        end
    end
end
Cons7 = [ {gamma1(:) >= 0 } ];
Cons8 = [ {gamma2(:) >= 0 } ];
Cons9 = [ {gamma3(:) >= 0 } ];
Cons10 = [ {(gamma1( find(ZeroGamma(:,:,:,:,1)==1) ) == 0)}] ;
Cons11 = [ {(gamma2( find(ZeroGamma(:,:,:,:,2)==1) ) == 0)}] ;
Cons12 = [ {(gamma3( find(ZeroGamma(:,:,:,:,3)==1) ) == 0)}] ;
Cons13 = [ {sum(x,2) + x_r <= 4*Capacity'} ];

Cons =  [Cons1;
        Cons2;
        Cons3;
        Cons4;
        Cons5;
        Cons6;
        Cons7;
        Cons8;
        Cons9;
        Cons10;
        Cons11;
        Cons12;
        Cons13;
        ];

%options = sdpsettings('verbose',1,'solver','Gurobi', 'showprogress',1);
opti.minimize(Obj);
opti.subject_to(Cons);
solver_options.expand = 1;
solver_options.print_time = 1;

%opti.solver('ipopt', solver_options);
opt = opti.save();
%Opti_RP = optimizer(Cons,Obj,options,{Demm},{Obj,gamma,x,x_r,D,epsilon});

% 
% for ii=1:N_nodes
%     Demands(ii,ii) = -sum(Demands(:,ii)) - Demands(ii,ii);
% end
% %%
% 
% sol1 = Opti_RP(Demands);
% %sol = optimize(Cons,Obj,options);
% sol.x = sol1{3};
% sol.xr = sol1{4};
% sol.obj = sol1{1};
% %sol.time =  G_road.Edges.Weight';
% %x_mat = sol.x;
% %x_mat(x_mat~=0) = 1;
% %sol.IndividualTimes = G_road.Edges.Weight'*x_mat;
% sol.Dem = sol1{5};
% sol.Gamma = sol1{2};
% sol.avgtime = sol.obj/TotDems;
% sol.eps = sol1{6};
% sol.Flows = sum(sol.x,2) + sol.xr;
% %sol.usertime = G_road.Edges.Weight'* sum( value(x) ,2);
% %sol.rebtime = G_road.Edges.Weight'*value(x_r);
% end

