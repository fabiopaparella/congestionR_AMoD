function [sol] = LTIFM_reb_cong(Demands,city)
load(strcat(city,'/Graphs.mat'));
Adj = adjacency(G_road);
Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);

for ii=1:N_nodes
    Demands(ii,ii) = -sum(Demands(:,ii)) - Demands(ii,ii);
end

x = sdpvar(N_edges*N_nodes,1,'full');
x_r = sdpvar(N_edges,1,'full');

B_kron = kron(eye(N_nodes),Binc);
FFT=kron(ones(1,N_nodes),G_road.Edges.Weight'); 
b=reshape(Demands,[],1);
Capacity = 300*ones(1,N_edges);
tau = G_road.Edges.Weight'./ Capacity; 
Capacity_th = 1.2*Capacity;
epsilon  = max(0,(sum( reshape(x,N_edges,[]) , 2) + x_r - Capacity_th' ));


Obj =  G_road.Edges.Weight'.*( sum( reshape(x,N_edges,[]) , 2) + tau* ( epsilon.^2 + epsilon.*(Capacity_th - x_r))); 
%Obj = G_road.Edges.Weight'.*(1+ (sum( reshape(x,N_edges,[]) , 2) + x_r)./Capacity' ).^4*(sum( reshape(x,N_edges,[]) , 2) + x_r);
Cons = [];
Cons1 = [ B_kron*x == b ];
Cons2 = [ x >=  0  ];
          %x <= 10000 ];
Cons3 = [ x_r >= 0 ];
          %x_r <= 10000];
      
Cons4 = [ Binc * ( (sum ( reshape(x,N_edges,[]) ,2) + x_r))  == 0   ];

Cons = [Cons1
        Cons2
        Cons3
        Cons4 ];

options = sdpsettings('verbose',1,'solver','gurobi', 'showprogress',1);
sol1 = optimize(Cons,Obj,options);
sol.x=value(x);
sol.xr =value(x_r);
sol.obj = FFT*value(x) + G_road.Edges.Weight'*value(x_r);
sol.time = G_road.Edges.Weight';
x_mat = reshape(sol.x,N_edges,N_nodes);
x_mat(x_mat~=0) = 1;
sol.IndividualTimes = G_road.Edges.Weight'*x_mat;
sol.Dem = Demands;
sol.usertime = FFT*value(x);
sol.rebtime = G_road.Edges.Weight'*value(x_r);
end

