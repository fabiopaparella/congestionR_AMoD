function [sol1] = TAP_casadi2(Demands,city,ExogT)
import casadi.*
opti = casadi.Opti()

load(strcat(city,'/Graphs.mat'));
Adj = adjacency(G_road);
Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);
TotDems = sum(Demands,'all');
for ii=1:N_nodes
    Demands(ii,ii) = -sum(Demands(:,ii)) - Demands(ii,ii);
end
if strcmp('SF', city)
    Capacity = Capacity/24;
else
    Capacity = 2000*ones(1,N_edges);
end
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


%% while loop
opti.set_value(ExogT_par, ExogT);


sol = opti.solve();
sol1.x = opti.value(x);
sol1.Flows = opti.value(y);
sol1.obj = opti.value(Obj);

end
