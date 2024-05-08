%% Code for experiment C. Collabrative agents in the paper Simplicial vector autoregressive models
addpath SC_VAR/
addpath Dependencies/
addpath Benchmark_algorithms/
addpath Football_data/
% Tstep=1 and Tstep=6 are used in the paper.

clc; clear all


T=500; % Total length of time series
P=6;   % P time lagged values are used to predict future
Tstep=6; % use t-P:t-1 data to predict time series vaule at time t+6

%% Intialize SC-VAR, S-VAR
param.T=T; % Total length of time series
param.Tstep=Tstep; % use t-P:t-1 data to predict time series vaule at time t+6
param.P=P;   % P time lagged values are used to predict future
scale=4; % scaling paramter for the simplicial convolution filters
param.K0=scale*1; % filter order L_0 filter
param.K1=scale*2 ;% filter order L_1 filter
param.K2=scale*1;% filter order L_2 filter
param.mu_0=0.01; % Regularization parameter for L_0
param.mu_1l=0.0000; %Regularization parameter for L1_lower  
param.mu_1u=0.0;  %Regularization parameter for L1_upper
param.mu_2=0.0; %Regularization parameter for L2
param.lambda=.01; % Regularization parameter for Beta
param.LassoEn=0; % LASSO regularizer enabler 
param.IncidDisplay=0;% Display incidence matrix
param.DisplayGraph=1; %Display the net work
param.HodgeNormlzn=0; %Enable Hodge laplacian normazlization 
param.FeatureNormlzn=1; % Enable feature normalization 
param.b=1; % normalization parameter
param.gamma=0.98; % Weight of forgetting window in RLS
param.delta=1-param.gamma;



param.edg=[3,10;10,6;6,4;3,2;6,1;6,7;4,7;4,5;5,7;7,1;1,2;2,9;9,8;8,5]


param.edg_index=1:size(param.edg,1);

Hodge=Topology_Generator(param);





for k=1:1 %loop to repeat experiment 
     signal_node=csvread('FB_signal_node.csv');% Read node data, speed of players in this case
     signal_edge=csvread('FB_signal_edge.csv'); % Read edge data,  distance between players
     signal_tri=csvread('FB_signal_tri.csv');% Read edge data,  area between players
     
   [n_node,nN]=size(signal_node); %extract number of nodes 
   [n_edge,nE]=size(signal_edge);%extract number of edges
   [n_tri,nT]=size(signal_tri); %extract number of triangles
   signal_compact=[signal_node(:,1:T);signal_edge(:,1:T);signal_tri(:,1:T)]; % stack the signals for other algorithms


%% SC-VAR 
param.V_En=1; % Vertex information enabler
param.F_En=1;  % Edge information enabler
param.T_En=1;  % Triangle information enabler
%call the function SC_VAR to get prediction NMSE :nmse_node_SC
[nmse_node_SC(k,:)]=SC_VAR(signal_node(:,1:T),signal_edge(:,1:T),signal_tri(:,1:T),Hodge,param);




%% S-VAR
param.V_En=1; % Vertex information enabler
param.F_En=0; % Edge information enabler
param.T_En=0; % Triangle information enabler
%call the function S_VAR to get prediction NMSE :nmse_node_S
[nmse_node_S(k,:)]=SC_VAR(signal_node(:,1:T),signal_edge(:,1:T),signal_tri(:,1:T),Hodge,param);

%% Beanch mark algorithms to for comparison
%%%%%%%%%%%%%% TIRSO %%%%%%%%%%%%%%%
[n_row_mx,nTimeInstants]=size(signal_compact);
tirsoObj = Tirso; % set tirso object up
tirsoObj.noOfNodes = n_row_mx;
tirsoObj.order     = P; 
tirsoObj.regPar    = 0;
tirsoObj.b_shrinkSelfLoops  = 0; % Bolstad
tirsoObj.forgettingFactor   = 0.98;
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
 nmse_node_tirso(k,:)=mean(CompNMSE(signal_node(:,1:T),m_prediction(1:n_node,1:end-Tstep)));
 nmse_edge_tirso(k,:)=mean(CompNMSE(signal_edge(:,1:T),m_prediction(n_node+1:n_node+n_edge,1:end-Tstep)));
 nmse_triangle_tirso(k,:)=mean(CompNMSE(signal_tri(:,1:T),m_prediction(n_node+n_edge+1:end,1:end-Tstep)));

 %%%%%%%%%%%%%% RFNL-TIRSO %%%%%%%%%%%%%%%

 

RFObj = RF_nltirso; % set tirso object up
RFObj.noOfNodes = n_row_mx;
RFObj.filtOrder = 6; % 
RFObj.lambda    = 1/10;
RFObj.NoOfRF    =10;
RFObj.vsigma    =.5*ones(RFObj.noOfNodes,1);

RFObj.forgettingFactor=.98 ;
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
  nmse_node_rfnltirso(k,:)=mean(CompNMSE(signal_node(:,1:T),RF_m_prediction(1:n_node,1:end-Tstep)));
 nmse_edge_rfnltirso(k,:)=mean(CompNMSE(signal_edge(:,1:T),RF_m_prediction(n_node+1:n_node+n_edge,1:end-Tstep)));
 nmse_triangle_rfnltirso(k,:)=mean(CompNMSE(signal_tri(:,1:T),RF_m_prediction(n_node+n_edge+1:end,1:end-Tstep)));
end

%% Plot the results
figure
subplot(3,1,1)
plot(signal_node');title('node signal' ) % plot node signal
subplot(3,1,2)
plot(signal_edge');title('edge signal') % plot edge signal
subplot(3,1,3)
plot(signal_tri');title('triangle signal')  % plot trangle signal

if param.V_En

figure
plot(mean(nmse_node_SC,1),'black-','MarkerIndices',1:20:500,'LineWidth',3,'DisplayName',strcat('SC-VAR :',num2str(Tstep),' step')) 

% Plot SC-VAR estimate in black color
hold on
plot(mean(nmse_node_tirso,1)','b-.','MarkerIndices',1:20:500,'LineWidth',3,'DisplayName',strcat('SC-VAR :',num2str(Tstep),' step')) % Plot TIRSO estimate in blue color
plot(mean(nmse_node_rfnltirso,1)','r.','MarkerIndices',1:5:500,'LineWidth',3,'DisplayName',strcat('SC-VAR :',num2str(Tstep),' step')) % Plot RFNL-TIRSO estimate in red color
plot(mean(nmse_node_S,1)','c--','MarkerIndices',1:20:500,'LineWidth',3,'DisplayName',strcat('SC-VAR :',num2str(Tstep),' step')) % Plot S-VAR estimate in cyan color
legend('show','FontSize', 18)
%legend("SC-VAR ","TIRSO","RFNL-TIRSO","S-VAR",'FontSize', 18)
%legend(sprintf('SC-VAR : %g step, TIRSO : %g step, RFNL-TIRSO : %g step, S-VAR : %g step', Tstep, Tstep,Tstep, Tstep),'FontSize', 18,'NumColumns',2);
xlabel('t','FontSize', 18)
ylabel('NMSE','FontSize', 18)
grid on


end


