city='SF';
Delay = 10;
load(strcat(city,'/Graphs.mat'));
Adj = adjacency(G_road);
Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);
set(gca,'ticklabelinterpreter','Latex','fontsize',16)
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultaxesticklabelinterpreter','latex')
set(groot,'defaultlegendinterpreter','latex')
vec = [2 5 10 15]; %[1 5 10 15]; 
Mar = ['o','+','*','x','v']; %'o','*','x','v'
colors=lines;
CM=[colors(1,:);colors(180,:);colors(220,:);colors(30,:);colors(256,:);]; %
A  = []; A_eps = []; AA=[];
B = []; B_eps = []; BB=[];
C = []; C_eps = []; CC=[];
multA = [1.5];
% for multiplicator = multA
%     A1=load(strcat(city,'/Congestion/Results/Mixed/solBase_Delay10Exo1','Dem',num2str(multiplicator),'.mat'));
%     A2 =load(strcat(city,'/Congestion/Results/Mixed/solBase_Delay10Exo0.5','Dem',num2str(multiplicator),'.mat'));
% 
%     A = [A A1.solBase.obj];
%     A_eps = [A_eps 4*A1.solBase.eps];%/max(A1.solBase.eps)];
%     AA = [AA A2.solBase.obj];
% end
%%
multB = [2];
for multiplicator = multB
    A1=load(strcat(city,'/Congestion/Results/solKn_Delay10Exo0.7','Dem',num2str(multiplicator),'.mat'));
    B1=load(strcat(city,'/Congestion/Results/Mixed/solKn_Delay10Exo0.7','Dem',num2str(multiplicator),'.mat'));
    C1=load(strcat(city,'/Congestion/Results/Mixed/solRP_Delay10Exo0.7','Dem',num2str(multiplicator),'.mat'));
   
    A_eps = [A_eps max(0,(A1.solKN.Flows + A1.solKN.Private - Capacity' )./Capacity')];
    B_eps = [B_eps B1.solKN.eps./Capacity'];
    C_eps = [C_eps C1.solRP.eps./Capacity'];
    
end
%%
% figure()
% hold on; grid on; box on;
% plot(multA,A./BB,'b')
% plot(multB,B./BB,'r')
% plot(multB,C./BB,'g')
% plot(multB,AA./BB,'k')
% plot(multB,CC./BB,'y')
%%
% figure()
% plot(multB,C./B)
%%
% figure()
% G_1 = G_road;
% A_eps_temp = A_eps(:,end);
% A_eps_temp(A_eps_temp==0)=0.001;
% G_1.Edges.Weight = A_eps_temp/1000;
% A_und = adjacency(G_1,'weighted');
% G_und = graph((A_und + A_und')/2);
% pp= plot(graph(Adj),'XData',NodesLoc(:,1),'YData',NodesLoc(:,2),'LineWidth',3)
% %G_und.Edges.EdgeColors = G_und.Edges.Weight;%A_eps(:,end);
% pp.EdgeCData = G_und.Edges.Weight;
% colormap(flipud(autumn))
% cc=colorbar('ticklabelinterpreter','Latex','fontsize',16)
% caxis([0,2.5]);
% cc.Label.String="Congestion";
% cc.Label.Interpreter = 'Latex';
% title('Ride-Sharing, $$\phi = 1$$')
% grid on; box on; 
% set(gca,'ticklabelinterpreter','Latex','fontsize',16,'XTick',[], 'YTick', [])

%%
% figure()
% G_1 = G_road;
% B_eps_temp = B_eps(:,1);
% B_eps_temp(B_eps_temp==0)=0.000001;
% G_1.Edges.Weight = B_eps_temp./Capacity';
% A_und = adjacency(G_1,'weighted');
% G_und = graph((A_und + A_und')/2);
% pp= plot(graph(Adj),'XData',NodesLoc(:,1),'YData',NodesLoc(:,2),'LineWidth',3)
% %G_und.Edges.EdgeColors = G_und.Edges.Weight;%A_eps(:,end);
% pp.NodeLabel = [];
% pp.EdgeCData = G_und.Edges.Weight;
% colormap(flipud(copper))
% cc=colorbar('ticklabelinterpreter','Latex','fontsize',16)
% caxis([0,3]);
% cc.Label.String="Congestion";
% cc.Label.Interpreter = 'Latex';
% %title('Unaware Assignments, $$\phi = 0.8$$')
% grid on; box on; 
% set(gca,'ticklabelinterpreter','Latex','fontsize',16,'XTick',[], 'YTick', [])

%%
figure()
G_1 = G_road;
C_eps(C_eps==0)=0.00000001;
G_1.Edges.Weight = C_eps;
A_und = adjacency(G_1,'weighted');
G_und = graph((A_und + A_und')/2);
pp= plot(graph(Adj),'XData',NodesLoc(:,1),'YData',NodesLoc(:,2),'LineWidth',3)
pp.NodeLabel = []
pp.EdgeCData = G_und.Edges.Weight;
colormap(flipud(copper))
cc=colorbar('ticklabelinterpreter','Latex','fontsize',16)
caxis([0,2.5])
cc.Label.String="Congestion";
cc.Label.Interpreter = 'Latex';
%title('Aware Assignments, $$\phi = 0.8$$')
grid on; box on; 
set(gca,'ticklabelinterpreter','Latex','fontsize',16,'XTick',[], 'YTick', [])

%%
subplot(1,2,1)
boxchart(B_eps-C_eps,'BoxFaceColor', CM(1,:))
set(gca,'ticklabelinterpreter','Latex','fontsize',16,'XTick',[])
grid on; box on;
ylim([-0.02,0.025])
%xlabel('Link Number')
ylabel('Congestion Difference')
set(gca,'ticklabelinterpreter','Latex','fontsize',16)
subplot(1,2,2)
boxchart(A_eps-C_eps,'BoxFaceColor', CM(1,:))
grid on; box on;
ylim([-0.5,1.7])
%yticklabels({''})
set(gca,'ticklabelinterpreter','Latex','fontsize',16,'XTick',[])

%%
figure()
subplot(1,2,1)
plot(A_eps-C_eps,'.','MarkerSize',15,'Color', CM(1,:))
grid on; box on;
ylim([-0.5,2])
xlabel('Link Number')
%ylabel('$$\kappa$$ Unaware - $$\kappa$$ Aware')
set(gca,'ticklabelinterpreter','Latex','fontsize',16)
subplot(1,2,2)
boxchart(A_eps-C_eps,'BoxFaceColor', CM(1,:))
grid on; box on;
ylim([-0.5,2])
%yticklabels({''})
set(gca,'ticklabelinterpreter','Latex','fontsize',16,'XTick',[])
