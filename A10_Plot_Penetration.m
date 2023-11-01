clc
clear all
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
Mar = ['o','+','*','x','v']; %'o','*','x','v'
mult = [1 1.5]; %
colors=lines;
CM=[colors(1,:);colors(180,:);colors(220,:);colors(30,:);colors(256,:);]; %
A  = []; A_eps = []; A_rs=[]; A_priv = []; D = []; E = [];
B = []; B_eps = []; B_rs=[]; B_priv = [];  
C = []; C_eps = []; C_rs=[]; C_priv = [];  
% for PenRate = [0.01 0.1:0.1:0.9 0.99] 
%     A_rs =[];
%     A_priv =[];
%    for multiplicator = mult
%     A1 = load(strcat(city,'/Congestion/Results/Mixed/solBase_Delay10Exo',num2str(PenRate),'Dem',num2str(multiplicator),'.mat'));
%     Flow = A1.solBase.Flows + A1.solBase.Private;
%     TravelTime = TRB(Flow,city);
%     A_rs = [A_rs TravelTime*sum(A1.solBase.x,2)/(multiplicator*PenRate*15000)];
%     A_priv = [A_priv TravelTime*A1.solBase.Private/(multiplicator*(1-PenRate)*15000)];
%    end
%  A = [A; A_rs A_priv];
%  
% end
for PenRate = [0.5:0.1:1]
    B_rs =[]; C_rs =[]; B_priv = []; C_priv = [];
for multiplicator = mult
    B1 =load(strcat(city,'/Congestion/Results/Mixed/solKn_Delay10Exo',num2str(PenRate),'Dem',num2str(multiplicator),'.mat'));
    C1 =load(strcat(city,'/Congestion/Results/Mixed/solRP_Delay10Exo',num2str(PenRate),'Dem',num2str(multiplicator),'.mat'));
    Flow1 = B1.solKN.Flows + B1.solKN.Private;
    Flow2 = C1.solRP.Flows + C1.solRP.Private;
    TravelTime1 = TRB(Flow1,city);
    TravelTime2 = TRB(Flow2,city);    
    B_rs = [B_rs (2*B1.solKN.RatioRP + 1 - B1.solKN.RatioRP)*TravelTime1*sum(B1.solKN.x,2)/(multiplicator*PenRate*15000)];
    B_priv = [B_priv TravelTime1*B1.solKN.Private/(multiplicator*(1-PenRate)*15000)];
    C_rs = [C_rs 2*TravelTime2*C1.solRP.Flows/(multiplicator*PenRate*15000) ];
    C_priv = [C_priv TravelTime2*C1.solRP.Private/(multiplicator*(1-PenRate)*15000)];
    %B_rs = [B_rs B1.solKN.obj/(multiplicator*PenRate*15000)];
    %B_priv = [B_priv TravelTime1*B1.solKN.Private/(multiplicator*(1-PenRate)*15000)];
    %C_rs = [C_rs C1.solRP.obj/(multiplicator*PenRate*15000)];
    %C_priv = [C_priv TravelTime2*C1.solRP.Private/(multiplicator*PenRate*15000)];
   
end
B = [B; B_rs]% B_priv]; 
C = [C; C_rs]% C_priv];
D = [D; B_rs];
E = [E; C_priv];
end
C(C<0.01) = inf;
%%
PenRate = [0.5:0.1:1];
figure()
hold on; grid on; box on;
plot(PenRate,B(:,1),'Color',CM(1,:),'Marker','o','LineWidth',1.5)
plot(PenRate,B(:,2),'Color',CM(2,:),'Marker','o','LineWidth',1.5)
plot(PenRate,C(:,1)','Color',CM(3,:),'Marker','x','LineWidth',1.5)
plot(PenRate,C(:,2),'Color',CM(4,:),'Marker','x','LineWidth',1.5)

plot(PenRate,D(:,1),'Color',CM(1,:),'Marker','o','LineWidth',1.5)
plot(PenRate,D(:,2),'Color',CM(2,:),'Marker','o','LineWidth',1.5)
plot(PenRate,E(:,1)','Color',CM(3,:),'Marker','x','LineWidth',1.5)
plot(PenRate,E(:,2),'Color',CM(4,:),'Marker','x','LineWidth',1.5)
% h(1)=plot(NaN, NaN,'Color','r','LineWidth',2);
% h(2)=plot(NaN, NaN,'Color','b','LineWidth',2);
% h(3)=plot(NaN, NaN,'Color','k','LineWidth',2);
%title('15000 Demands per Hour ')
% legend('AMoD (RS)','Private (RS)','AMoD (RP)','Private (RP)')%,'7500 Demands/h','15000 Demands/h','NumColumns',2)
legend('Unaware RP 15k Requests/h','Unaware RP 22.5k Requests/h','Aware RP 15k Requests/h','Aware RP 22.5k Requests/h')%,'7500 Demands/h','15000 Demands/h','NumColumns',2)
set(gca,'ticklabelinterpreter','Latex','fontsize',14)
ylabel('Average Travel Time [min]')
xlabel('Penetration Rate $$\phi$$')


% figure()
% hold on; grid on; box on;
% plot(PenRate,A(:,2),'Color',CM(1,:),'Marker','o','LineWidth',1.5)
% plot(PenRate,A(:,4),'Color',CM(2,:),'Marker','o','LineWidth',1.5)
% 
% plot(PenRate,B(:,2),'Color',CM(3,:),'Marker','x','LineWidth',1.5)
% plot(PenRate,B(:,4),'Color',CM(2,:),'Marker','x','LineWidth',1.5)
% xlim([0,1])
% ylim([5 35])
% h(1)=plot(NaN, NaN,'Color','r','LineWidth',2);
% h(2)=plot(NaN, NaN,'Color','b','LineWidth',2);
% h(3)=plot(NaN, NaN,'Color','k','LineWidth',2);
% title('15000 Demands per Hour ')
% 
% set(gca,'ticklabelinterpreter','Latex','fontsize',14)
% ylabel('Average Travel Time [min]')
% xlabel('Penetration Rate $$\phi$$')
