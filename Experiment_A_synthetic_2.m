%% This code runs Experiment A. compariosn between SC-VAR and S-VAR for sythetic data
clc; close all; clear all


T=1000; % Total length of time series
P=4; %  P time lagged values are used to predict future
Tstep=1; % use t-P:t-1 data to predict time series vaule at time t


%% Intialize SC-VAR, S-VAR
param.Tstep=Tstep; % use t-P:t-1 data to predict time series vaule at time t+Tstep-1
param.T=T; % Total length of time series
param.P=P;% P time lagged values are used to predict future
scale=4 ;  % scaling paramter for the simplicial convolution filters
param.K0=scale*1; % filter order L_0 filter
param.K1=scale*2;  % filter order L_1 filter
param.K2=scale*1; % filter order L_2 filter
param.mu_0=0; % Regularization parameter for L_0
param.mu_1l=0.0000; %Regularization parameter for L1_lower
param.mu_1u=0.0; %Regularization parameter for L1_upper
param.mu_2=0.0;  %Regularization parameter for L2
param.lambda=.001; % Regularization parameter for Beta
param.LassoEn=0; % LASSO regularizer enabler 
param.IncidDisplay=0; % Display incidence matrix 0 or 1
param.DisplayGraph=0; %Display the net work 0 or 1
param.HodgeNormlzn=0; %Enable Hodge laplacian normazlization 
param.FeatureNormlzn=1;  % Enable feature normalization 
param.b=1; % normalization parameter
param.gamma=0.99; % Weight of forgetting window in RLS
param.delta=1-param.gamma;



 param.edg=[1,2;2,3;3,4;4,5;5,6;6,7;7,8;8,9;9,10;10,11;11,12;12,1;...
            2,4;4,6;6,8;8,10;10,12;12,2;...
            2,13;13,4;6,13;13,8;10,13;13,12];

param.edg_index=1:size(param.edg,1);

Hodge=Topology_Generator(param); %Generate SC structre





for k=1:5 % loop to repeat experiment
display(strcat("Experiment Num:",num2str(k)))
% function Simplicial_Signal_Generator return synthtic data for node, edge and triangle signal
   [signal_node,signal_edge,signal_tri]=Simplicial_Signal_Generator_2(Hodge,param);
   [n_node,nT]=size(signal_node);
   [n_edge,nT]=size(signal_edge);
   [n_tri,nT]=size(signal_edge);
   signal_compact=[signal_node;signal_edge;signal_tri];



param.V_En=1;
param.F_En=1;
param.T_En=1;
[nmse_node_SC(k,:),nmse_edge_SC(k,:),nmse_tri_SC(k,:),maxeig_V,maxeig_F,maxeig_T]=SC_VAR(signal_node,signal_edge,signal_tri,Hodge,param);

param.V_En=0;
param.F_En=1;
param.T_En=0;
[nmse_node_S_0(k,:),nmse_edge_S(k,:),nmse_tri_S_0(k,:)]=SC_VAR(signal_node,signal_edge,signal_tri,Hodge,param);




param.V_En=1;
param.F_En=0;
param.T_En=0;
[nmse_node_S(k,:),nmse_edge_S_1(k,:),nmse_tri_S_1(k,:)]=SC_VAR(signal_node,signal_edge,signal_tri,Hodge,param);

param.V_En=0;
param.F_En=0;
param.T_En=1;
[nmse_node_S_2(k,:),nmse_edge_S_2(k,:),nmse_tri_S(k,:)]=SC_VAR(signal_node,signal_edge,signal_tri,Hodge,param);


end

figure
subplot(3,1,1)
plot(signal_node');title('node signal')
subplot(3,1,2)
plot(signal_edge');title('edge signal')
subplot(3,1,3)
plot(signal_tri');title('triangle signal')
figure
plot(mean(nmse_node_SC,1),'black-*','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('SC-VAR :',num2str(Tstep),' step')) 
title('Vertex signal')
% Plot SC-VAR estimate in black color
hold on
plot(mean(nmse_node_S,1)','b-.*','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('S-VAR :',num2str(Tstep),' step')) % Plot S-VAR estimate in blue color
legend('show','FontSize', 18)
%legend("SC-VAR ","TIRSO","RFNL-TIRSO","S-VAR",'FontSize', 18)
%legend(sprintf('SC-VAR : %g step, TIRSO : %g step, RFNL-TIRSO : %g step, S-VAR : %g step', Tstep, Tstep,Tstep, Tstep),'FontSize', 18,'NumColumns',2);
xlabel('t','FontSize', 18)
ylabel('NMSE','FontSize', 18)
grid on

figure
plot(mean(nmse_edge_SC,1),'black-o','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('SC-VAR :',num2str(Tstep),' step')) 
title('Edge signal')
% Plot SC-VAR estimate in black color
hold on
plot(mean(nmse_edge_S,1)','b-.o','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('S-VAR :',num2str(Tstep),' step')) % Plot S-VAR estimate in blue color

legend('show','FontSize', 18)
%legend("SC-VAR ","TIRSO","RFNL-TIRSO","S-VAR",'FontSize', 18)
%legend(sprintf('SC-VAR : %g step, TIRSO : %g step, RFNL-TIRSO : %g step, S-VAR : %g step', Tstep, Tstep,Tstep, Tstep),'FontSize', 18,'NumColumns',2);
xlabel('t','FontSize', 18)
ylabel('NMSE','FontSize', 18)
grid on



figure
plot(mean(nmse_tri_SC,1),'black->','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('SC-VAR :',num2str(Tstep),' step')) 
title('Triangle signal')
% Plot SC-VAR estimate in black color
hold on
plot(mean(nmse_tri_S,1)','b-.>','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('S-VAR :',num2str(Tstep),' step')) % Plot S-VAR estimate in blue color
legend('show','FontSize', 18)
%legend("SC-VAR ","TIRSO","RFNL-TIRSO","S-VAR",'FontSize', 18)
%legend(sprintf('SC-VAR : %g step, TIRSO : %g step, RFNL-TIRSO : %g step, S-VAR : %g step', Tstep, Tstep,Tstep, Tstep),'FontSize', 18,'NumColumns',2);
xlabel('t','FontSize', 18)
ylabel('NMSE','FontSize', 18)
grid on


