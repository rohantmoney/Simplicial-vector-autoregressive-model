%% This code runs Experiment A. compariosn between SC-VAR aginest Benchmark algorithms for sythetic data
clc; close all; clear all
 addpath SC_VAR/
 addpath Dependencies/
 addpath Benchmark_algorithms/
 addpath Synthtic_data_generator/


T=2000; % Total length of time series
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





for k=1:25 % loop to repeat experiment
display(strcat("Experiment Num:",num2str(k)))
% function Simplicial_Signal_Generator return synthtic data for node, edge and triangle signal
   [signal_node,signal_edge,signal_tri]=Simplicial_Signal_Generator(Hodge,param);
   
   [n_node,nT]=size(signal_node);
   [n_edge,nT]=size(signal_edge);
   [n_tri,nT]=size(signal_edge);
   signal_compact=[signal_node;signal_edge;signal_tri];


if 1
param.V_En=1;
param.F_En=1;
param.T_En=1;


%SC_VAR function return the NMSE for triangle, edge and vertex signals
[nmse_node_SC(k,:),nmse_edge_SC(k,:),nmse_tri_SC(k,:)]=SC_VAR(signal_node,signal_edge,signal_tri,Hodge,param);





%%Benchmark algorithms to compare
%%%%%%%%%%%%%%%%%TIRSO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n_row_mx,nTimeInstants]=size(signal_compact);
tirsoObj = Tirso; % set tirso object up
tirsoObj.noOfNodes = n_row_mx;
tirsoObj.order     = P; % we can try a higher order later
tirsoObj.regPar    = 0;
tirsoObj.b_shrinkSelfLoops  = 0; % Bolstad
tirsoObj.forgettingFactor   = 0.99;
tirsoObj.h_stepsize         = @(ts)1/eigs(ts.m_Phi,1);
% initialize
tState_in = tirsoObj.initialize(0, signal_compact( :,1:tirsoObj.order)');
e_tirso=zeros(24,nTimeInstants);
tic
for t = tirsoObj.order+1:nTimeInstants
    mtemp= signal_compact(:, t);
    tState_in = tirsoObj.update(tState_in, mtemp);
    m_predic(:,:)=tState_in.predictManyFromBuffer(1)';
    m_prediction(1:tirsoObj.noOfNodes,t+1)= m_predic(:,1);
   
    
end
time_tirso=toc;
 nmse_node_tirso(k,:)=mean(CompNMSE(signal_node,m_prediction(1:n_node,1:end-1)));
 nmse_edge_tirso(k,:)=mean(CompNMSE(signal_edge,m_prediction(n_node+1:n_node+n_edge,1:end-1)));
 nmse_triangle_tirso(k,:)=mean(CompNMSE(signal_tri,m_prediction(n_node+n_edge+1:end,1:end-1)));
 %%%%%%%%%%%%%%%%RF_NLtiso%%%%%%%%%%%%%%%%


RFObj = RF_nltirso; % set tirso object up
RFObj.noOfNodes = n_row_mx;
RFObj.filtOrder = 4; % we can try a higher order later
RFObj.lambda    = 1/10000;
RFObj.NoOfRF    =7;
RFObj.vsigma    =.5*ones(RFObj.noOfNodes,1);

RFObj.forgettingFactor=.99 ;
RFObj.h_stepsize= @(RF_ts)1/eigs(RF_ts.m_Phi,1);
RFObj.eta       =400;
RFState_in = RFObj.initialize(1,signal_compact( :,1:RFObj.filtOrder)');
tic
 for t = RFObj.filtOrder+1:nTimeInstants
     
    mtemp= signal_compact(:, t);
    RFState_in = RFObj.update(RFState_in, mtemp);
    RF_m_predic(:,:)=RFState_in.predictManyFromBuffer(1)';
    RF_m_prediction(1:tirsoObj.noOfNodes,t+1)= RF_m_predic(:,1);
 end
 time_rf=toc;
  nmse_node_rfnltirso(k,:)=mean(CompNMSE(signal_node,RF_m_prediction(1:n_node,1:end-1)));
 nmse_edge_rfnltirso(k,:)=mean(CompNMSE(signal_edge,RF_m_prediction(n_node+1:n_node+n_edge,1:end-1)));
 nmse_triangle_rfnltirso(k,:)=mean(CompNMSE(signal_tri,RF_m_prediction(n_node+n_edge+1:end,1:end-1)));
end
%%%%%%%%%%%%%%%%%%moving average%%%%%%%%%%%%%%%%%%%
tic
for t = tirsoObj.order+1:nTimeInstants
    mv_prediction(1:tirsoObj.noOfNodes,t+1)=mean(signal_compact(:,1:t),2);
    
end
time_mav=toc;
 nmse_node_mav(k,:)=mean(CompNMSE(signal_node,mv_prediction(1:n_node,1:end-1)));
 nmse_edge_mav(k,:)=mean(CompNMSE(signal_edge,mv_prediction(n_node+1:n_node+n_edge,1:end-1)));
 nmse_triangle_mav(k,:)=mean(CompNMSE(signal_tri,mv_prediction(n_node+n_edge+1:end,1:end-1)));

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
plot(mean(nmse_node_tirso,1)','b-.*','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('TIRSO :',num2str(Tstep),' step')) % Plot TIRSO estimate in blue color
plot(mean(nmse_node_rfnltirso,1)','R:*','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('RFNL-TIRSO :',num2str(Tstep),' step')) % Plot RFNL-TIRSO estimate in red color
plot(mean(nmse_node_mav,1)','y*','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('Moving average :',num2str(Tstep),' step')) % Plot S-VAR estimate in yellow color
legend('show','FontSize', 18)

xlabel('t','FontSize', 18)
ylabel('NMSE','FontSize', 18)
grid on

figure
plot(mean(nmse_edge_SC,1),'black-o','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('SC-VAR :',num2str(Tstep),' step')) 
title('Edge signal')
% Plot SC-VAR estimate in black color
hold on
plot(mean(nmse_edge_tirso,1)','b-.o','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('TIRSO :',num2str(Tstep),' step')) % Plot TIRSO estimate in blue color
plot(mean(nmse_edge_rfnltirso,1)','r:o','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('RFNL-TIRSO :',num2str(Tstep),' step')) % Plot RFNL-TIRSO estimate in red color
plot(mean(nmse_edge_mav,1)','yo','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('Moving average :',num2str(Tstep),' step')) % Plot S-VAR estimate in yellow color
legend('show','FontSize', 18)

xlabel('t','FontSize', 18)
ylabel('NMSE','FontSize', 18)
grid on



figure
plot(mean(nmse_tri_SC,1),'black->','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('SC-VAR :',num2str(Tstep),' step')) 
title('Triangle signal')
% Plot SC-VAR estimate in black color
hold on
plot(mean(nmse_triangle_tirso,1)','b-.>','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('TIRSO :',num2str(Tstep),' step')) % Plot TIRSO estimate in blue color
plot(mean(nmse_triangle_rfnltirso,1)','R:>','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('RFNL-TIRSO :',num2str(Tstep),' step')) % Plot RFNL-TIRSO estimate in red color
plot(mean(nmse_triangle_mav,1)','y>','MarkerIndices',1:50:2000,'MarkerSize',15,'LineWidth',3,'DisplayName',strcat('Moving average :',num2str(Tstep),' step')) % Plot S-VAR estimate in yellow color
legend('show','FontSize', 18)
xlabel('t','FontSize', 18)
ylabel('NMSE','FontSize', 18)
grid on


