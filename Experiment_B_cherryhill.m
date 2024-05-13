%% Code for experiment B. Water distribution network in the paper Simplicial vector autoregressive models
% Tstep=1 and Tstep=2 are used in the paper.
clc; clear all
 addpath SC_VAR/
 addpath Dependencies/
 addpath Benchmark_algorithms/
 addpath Water_data/


T=112; % Total length of time series
P=5;   %  P time lagged values are used to predict future
Tstep=1; % use t-P:t-1 data to predict time series vaule at time t





%% Intialize SC-VAR, S-VAR
param.T=T; % Total length of time series
param.Tstep=Tstep; % use t-P:t-1 data to predict time series vaule at time t+Tstep-1
param.P=P;% P time lagged values are used to predict future
scale=2 ; % scaling paramter for the simplicial convolution filters
param.K0=scale*1; % filter order L_0 filter
param.K1=scale*2; % filter order L_1 filter
param.K2=scale*1;% filter order L_2 filter
param.mu_0=0.5; % Regularization parameter for L_0
param.mu_1l=.5000000; %Regularization parameter for L1_lower  
param.mu_1u=0.0; %Regularization parameter for L1_upper
param.mu_2=0.0; %Regularization parameter for L2
param.lambda=0.0001; % Regularization parameter for Beta
param.LassoEn=0.00001; % LASSO regularizer enabler 
param.IncidDisplay=0.0;% Display incidence matrix 0 or 1
param.DisplayGraph=0;  %Display the net work 0 or 1
param.HodgeNormlzn=0; %Enable Hodge laplacian normazlization 
param.FeatureNormlzn=1; % Enable feature normalization 
param.b=1;% normalization parameter

param.gamma=0.99; % Weight of forgetting window in RLS
param.delta=1-param.gamma;




directed_links=[1,2;2,5;2,3;3,4;4,5;5,6;6,7;7,8;7,9;8,10;9,11;11,12;12,13;13,14; 13,16;14,15;14,20;15,17;15,24;16,17;16,19;17,18;18,32;20,21;20,22;21,22;22,33;23,25;24,23; 25,26;25,31;27,29;28,35;28,36;29,28;29,35;31,27;32,19;33,34;35,30];% Edge set
[edg,sort_idx]=sortrows(directed_links); % Lexigocraphic sorting

 param.edg=edg;

param.edg_index=1:size(param.edg,1);

Hodge=Topology_Generator(param);  % generate the SC structure









for k=1:1  %loop to repeat experiment 
display(strcat("Experiment Num:",num2str(k)))
signal_node=csvread('Pressure.csv');% Read vertex data, pressure in this case
signal_edge=csvread('Flow.csv'); % Read edge data,  flow rate in this case
signal_tri=zeros(2,112); % traingle signals zeros

   
   [n_node,nT]=size(signal_node); %number of nodes
   [n_edge,nT]=size(signal_edge); %number of edges
   
   signal_compact=[signal_node;signal_edge; signal_tri];% stacked signal for other bench mark algorithms

%% SC-VAR
param.V_En=1; %Enable vertex data
param.F_En=1; %Enable edge data
param.T_En=0;
%call the function SC_VAR to get return NMSE :nmse_node_S
[nmse_node_SC(k,:),nmse_edge_SC(k,:)]=SC_VAR(signal_node,signal_edge,signal_tri,Hodge,param);

%% S-VAR
param.V_En=0;
param.F_En=1; %Enable edge data
param.T_En=0;
%call the function S_VAR to get return NMSE :nmse_node_S
[nmse_node_S(k,:),nmse_edge_S(k,:)]=SC_VAR(signal_node,signal_edge,signal_tri,Hodge,param);


%%
param.V_En=0;
param.F_En=1;
param.T_En=0;
tic
[nmse_node_010(k,:),nmse_edge_010(k,:),maxeig_V,maxeig_F]=SC_VAR(signal_node,signal_edge,signal_tri,Hodge,param);
tim_scvar=toc;
%% Beanch mark algorithms to for comparison

%%%%%%%%%%%%%%TIRSO%%%%%%%%%%%%%%%%
[n_row_mx,nTimeInstants]=size(signal_compact);
tirsoObj = Tirso; % set tirso object up
tirsoObj.noOfNodes = n_row_mx;
tirsoObj.order     = P; % we can try a higher order later
tirsoObj.regPar    = 0.000008;
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
    m_predic(:,:)=tState_in.predictManyFromBuffer(Tstep)';
    m_prediction(1:tirsoObj.noOfNodes,t+Tstep)= m_predic(:,Tstep);
   
    
end
time_tirso=toc;
 nmse_node_tirso(k,:)=mean(CompNMSE(signal_node,m_prediction(1:n_node,1:end-Tstep)));
 nmse_edge_tirso(k,:)=mean(CompNMSE(signal_edge,m_prediction(n_node+1:n_node+n_edge,1:end-Tstep)));
%%%%%%%%%%%%%%RFNL-TIRSO%%%%%%%%%%%%%%%%


RFObj = RF_nltirso; % set tirso object up
RFObj.noOfNodes = n_row_mx;
RFObj.filtOrder = 4; % we can try a higher order later
RFObj.lambda    = 1/1000;
RFObj.NoOfRF    =5;
RFObj.vsigma    =1*ones(RFObj.noOfNodes,1);

RFObj.forgettingFactor=.99 ;
RFObj.h_stepsize= @(RF_ts)1/eigs(RF_ts.m_Phi,1);
RFObj.eta       =200;
RFState_in = RFObj.initialize(1,signal_compact( :,1:RFObj.filtOrder)');
tic
 for t = RFObj.filtOrder+1:nTimeInstants
     
    mtemp= signal_compact(:, t);
    RFState_in = RFObj.update(RFState_in, mtemp);
    RF_m_predic(:,:)=RFState_in.predictManyFromBuffer(Tstep)';
    RF_m_prediction(1:tirsoObj.noOfNodes,t+Tstep)= RF_m_predic(:,Tstep);

      
 end
 time_rf=toc;
  nmse_node_rfnltirso(k,:)=mean(CompNMSE(signal_node,RF_m_prediction(1:n_node,1:end-Tstep)));
 nmse_edge_rfnltirso(k,:)=mean(CompNMSE(signal_edge,RF_m_prediction(n_node+1:n_node+n_edge,1:end-Tstep)));
end

%% Plot the results
figure
subplot(3,1,1)
plot(signal_node');title('node signal')
subplot(3,1,2)
plot(signal_edge');title('edge signal')


figure
plot(mean(nmse_edge_SC,1),'black-','MarkerIndices',1:10:100,'LineWidth',3,'DisplayName',strcat('SC-VAR :',num2str(Tstep),' step')) 

% Plot SC-VAR estimate in black color
hold on
plot(mean(nmse_edge_tirso,1)','b-.','MarkerIndices',1:10:110,'LineWidth',3,'DisplayName',strcat('TIRSO :',num2str(Tstep),' step')) % Plot TIRSO estimate in blue color
plot(mean(nmse_edge_rfnltirso,1)','r.','MarkerIndices',1:2:110,'LineWidth',3,'DisplayName',strcat('RFNL-TIRSO :',num2str(Tstep),' step')) % Plot RFNL-TIRSO estimate in red color
plot(mean(nmse_edge_S,1)','c--','MarkerIndices',1:10:100,'LineWidth',3,'DisplayName',strcat('S-VAR :',num2str(Tstep),' step')) % Plot S-VAR estimate in cyan color
legend('show','FontSize', 18)
%legend("SC-VAR ","TIRSO","RFNL-TIRSO","S-VAR",'FontSize', 18)
%legend(sprintf('SC-VAR : %g step, TIRSO : %g step, RFNL-TIRSO : %g step, S-VAR : %g step', Tstep, Tstep,Tstep, Tstep),'FontSize', 18,'NumColumns',2);
xlabel('t','FontSize', 18)
ylabel('NMSE','FontSize', 18)
grid on


