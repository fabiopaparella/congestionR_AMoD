clc
clear all
city='NYC120';
Delay = 10;
load(strcat(city,'/Graphs.mat'));
Adj = adjacency(G_road);
Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);
set(gca,'ticklabelinterpreter','Latex','fontsize',16)
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultaxesticklabelinterpreter','latex')
set(groot,'defaultlegendinterpreter','latex')
mult = [4]; %
colors=lines;
CM=[colors(1,:);colors(180,:);colors(220,:);colors(30,:);colors(256,:);]; %
A  = []; A_eps = []; A_rs=[]; A_priv = []; A_rp = []; EPS = [];
B = []; B_eps = []; B_rs=[]; B_priv = [];  B_rp = [];
C = []; C_eps = []; C_rs=[]; C_priv = [];  C_rp = [];
for PenRate = [0.01 0.1:0.1:0.9 0.99]
    A_rs =[]; A_rp = []; A_priv =[]; EPs = [];
    for Psi = [0.01 0.1:0.1:0.9 0.99]
        multiplicator = mult;
        A1 = load(strcat(city,'/Congestion/Results/3Mixed/Delay10Exo',num2str(PenRate),'RP',num2str(Psi),'Dem',num2str(multiplicator),'.mat'));
        Flow = A1.solBase.Flows + A1.solKN.Flows + A1.privateBase.Flows;
        TravelTime = TRB(Flow,city);
        A_rs = [A_rs TravelTime*sum(A1.solBase.x,2)/(multiplicator*PenRate*(1-Psi)*53000)];
        A_rp = [A_rp (2*A1.solKN.RatioRP + 1 - A1.solKN.RatioRP)*TravelTime*sum(A1.solKN.x,2)/(multiplicator*PenRate*Psi*53000)];
        A_priv = [A_priv TravelTime*A1.privateBase.Flows/(multiplicator*(1-PenRate)*53000)];
        EPs = [EPs mean(A1.solBase.eps)];
    end
    A = [A; A_rs];
    B = [B ; A_rp];
    C = [C; A_priv];
    EPS = [EPS; EPs];
end

%%
Phi = [0:0.1:1];
psi = [0:0.1:1];
[XX,YY]=meshgrid(Phi,psi);
fig=figure()
hold on; grid on; box on;
contourf(XX,YY,A')
caxis([10,18]);
%title('Avg. Travel Time RS Users [min] ')
set(gca,'ticklabelinterpreter','Latex','fontsize',18)
ylabel('Ride-pooling Penetration Rate $$\psi$$')
%xlabel('AMoD Penetration Rate $$\phi$$')
h = axes(fig,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on'; 
cc=colorbar(h,'ticklabelinterpreter','Latex','fontsize',18)
%cc.Label.String="Average Travel Time [min]";
cc.Position= [0.925 0.11 0.025 0.82];
cc.Label.Interpreter = 'Latex';
caxis([10,18]);
%%
fig=figure()
hold on; grid on; box on;
contourf(XX,YY,C')
caxis([10,18]);
%title('Avg. Travel Time Private Users [min]')
set(gca,'ticklabelinterpreter','Latex','fontsize',18)
%ylabel('Ride-pooling Penetration Rate $$\Psi$$')
%xlabel('AMoD Penetration Rate $$\phi$$')

h = axes(fig,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on'; 
cc=colorbar(h,'ticklabelinterpreter','Latex','fontsize',18)
%cc.Label.String="Average Travel Time [min]";
cc.Position= [0.925 0.11 0.025 0.82];
cc.Label.Interpreter = 'Latex';
caxis([10,18]);
%%
[XX,YY]=meshgrid(Phi,psi);
fig=figure()
grid on; box on;
contourf(XX,YY,B')
caxis([10,29]);
%title('Avg. Travel Time RP Users [min] ')
set(gca,'ticklabelinterpreter','Latex','fontsize',20)
ylabel('Ride-pooling Penetration Rate $$\psi$$')
xlabel('AMoD Penetration Rate $$\phi$$')
h = axes(fig,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on'; 
cc=colorbar(h,'ticklabelinterpreter','Latex','fontsize',18)
%cc.Label.String="Average Travel Time [min]";
cc.Position= [0.925 0.11 0.025 0.82];
cc.Label.Interpreter = 'Latex';
caxis([10,29]);
%%
[XX,YY]=meshgrid(Phi,psi);
fig=figure()

grid on; box on;
contourf(XX,YY,B'-A')
caxis([2,11]);
%title('Difference Between RP and RS [min]')
set(gca,'ticklabelinterpreter','Latex','fontsize',18)
%ylabel('Ride-pooling Penetration Rate $$\Psi$$')
xlabel('AMoD Penetration Rate $$\phi$$')

h = axes(fig,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on'; 
cc=colorbar(h,'ticklabelinterpreter','Latex','fontsize',18)
%cc.Label.String="Average Travel Time [min]";
cc.Position= [0.925 0.11 0.025 0.82];
cc.Label.Interpreter = 'Latex';
caxis([2,11]);

%%
fig=figure()
grid on; box on; hold on;

plot(NaN, NaN,'Color',CM(1,:),'LineWidth',2);
plot(NaN, NaN,'Color',CM(2,:),'LineWidth',2);
plot(NaN, NaN,'Marker','o','MarkerSize',8,'Color','k','Linestyle', 'none');
plot(NaN, NaN,'Marker','square','MarkerSize',8,'Color','k','Linestyle', 'none');
plot(NaN, NaN,'Marker','*','MarkerSize',8,'Color','k','Linestyle', 'none');

plot(Phi,A(:,3),'LineWidth',2,'Marker','o','MarkerSize',10,'Color',CM(1,:))
plot(Phi,B(:,3),'LineWidth',2,'Marker','square','MarkerSize',10,'Color',CM(1,:))
plot(Phi,C(:,3),'LineWidth',2,'Marker','*','MarkerSize',10,'Color',CM(1,:))

plot(Phi,A(:,9),'LineWidth',2,'Marker','o','MarkerSize',10,'Color',CM(2,:))
plot(Phi,B(:,9),'LineWidth',2,'Marker','square','MarkerSize',10,'Color',CM(2,:))
plot(Phi,C(:,9),'LineWidth',2,'Marker','*','MarkerSize',10,'Color',CM(2,:))

set(gca,'ticklabelinterpreter','Latex','fontsize',14)
xlabel('Penetration Rate $$\phi$$')
ylabel('Average Travel Time [min]')
legend('$$\psi=0.2$$','$$\psi=0.8$$','Individual Ride-sharing','Ride-pooling','Private')
ylim([0 30])

%%
Phi = [0:0.1:1];
psi = [0:0.1:1];
fig=figure()
grid on; box on; hold on;

plot(NaN, NaN,'Color',CM(3,:),'LineWidth',2);
plot(NaN, NaN,'Color',CM(4,:),'LineWidth',2);
plot(NaN, NaN,'Marker','o','MarkerSize',8,'Color','k','Linestyle', 'none');
plot(NaN, NaN,'Marker','square','MarkerSize',8,'Color','k','Linestyle', 'none');
plot(NaN, NaN,'Marker','*','MarkerSize',8,'Color','k','Linestyle', 'none');

plot(Phi,A(3,:),'LineWidth',2,'Marker','o','MarkerSize',10,'Color',CM(3,:))
plot(Phi,B(3,:),'LineWidth',2,'Marker','square','MarkerSize',10,'Color',CM(3,:))
plot(Phi,C(3,:),'LineWidth',2,'Marker','*','MarkerSize',10,'Color',CM(3,:))

plot(Phi,A(9,:),'LineWidth',2,'Marker','o','MarkerSize',10,'Color',CM(4,:))
plot(Phi,B(9,:),'LineWidth',2,'Marker','square','MarkerSize',10,'Color',CM(4,:))
plot(Phi,C(9,:),'LineWidth',2,'Marker','*','MarkerSize',10,'Color',CM(4,:))

set(gca,'ticklabelinterpreter','Latex','fontsize',14)
xlabel('Ride-pooling Penetration Rate $$\psi$$')
ylabel('Average Travel Time [min]')
legend('$$\phi=0.2$$','$$\phi=0.8$$','Individual Ride-sharing','Ride-pooling','Private')
ylim([0 30])