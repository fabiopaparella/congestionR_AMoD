function [sol] = LTIFM_NRP_congestion(Demands,city,ExogT)
load(strcat(city,'/Graphs.mat'));
Adj = adjacency(G_road);
Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);
TotDems = sum(Demands,'all');
for ii=1:N_nodes
    Demands(ii,ii) = -sum(Demands(:,ii)) - Demands(ii,ii);
end

x = sdpvar(N_edges,N_nodes,'full');
x_r = sdpvar(N_edges,1,'full');
Exog_traffic = ExogT;%*Capacity'; %0*ones(1,N_edges)';
tau = G_road.Edges.Weight'./Capacity; 
Capacity_th = 1.2*Capacity;
epsilon  = max(0,sum(x,2) + x_r + Exog_traffic - Capacity_th' );


Obj =  G_road.Edges.Weight'*(sum(x,2) + x_r) + 7.3*tau* ( epsilon.^2 + epsilon.*(Capacity_th' - Exog_traffic ));  %FFT*x  5.27 3.75 for x 3
Cons = [];
Cons1 = [ Binc*x == Demands ];
Cons2 = [ x >=  0  ];
Cons3 = [ x_r >= 0 ];
Cons4 = [ Binc * ( (sum ( x,2) + x_r))  == 0   ];
Cons5 = [ sum(x,2) + x_r + Exog_traffic <= 3*Capacity' ];
Cons = [Cons1
        Cons2
        Cons3
        Cons4
        Cons5
       ];

options = sdpsettings('verbose',1,'solver','gurobi', 'showprogress',1);
sol1 = optimize(Cons,Obj,options);
sol.x=value(x);
sol.xr =value(x_r);
sol.obj = G_road.Edges.Weight'*(sum(value(x),2)+ + value(x_r)) + tau* ( value(epsilon).^2 + value(epsilon).*(Capacity_th' - Exog_traffic ));
sol.FFtime = G_road.Edges.Weight;
sol.Flows = sum(value(x),2) + value(x_r);
sol.eps = value(epsilon);
sol.Dem = Demands;
sol.Private = ExogT;
end

