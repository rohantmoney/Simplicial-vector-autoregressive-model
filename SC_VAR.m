
function [nmse_node,nmse_edge,nmse_tri,maxeig_V,maxeig_F,maxeig_T]=SC_VAR(signal_node,signal_edge,signal_tri,Hodge,param);
% The function estiame the parameters of the SC-VAR model in an online way.
% At each time instant forcast Tstep ahead prediction and find the corresponding Normalized mean squared error 
%% Inputs of function
% signal_vertex
% signal_edge
% signal_tri
% Hodge: conatin structure of network
% param: Hyper parameters and flags
%% Outputs of function
% NMSE: Tstep ahaed vertex signal
% NMSE: Tstep ahaed edge signal
% NMSE: Tstep ahaed triangle signal


%% Local copy of parameters 
T=param.T;  % Total length of time series
P=param.P; % P time lagged values are used to predict future
mu_0=param.mu_0; % Regularization parameter for L_0
mu_1l=param.mu_1l;  %Regularization parameter for L1_lower  
mu_1u=param.mu_1u;  %Regularization parameter for L1_upper
mu_2=param.mu_2; %Regularization parameter for L2
Tstep=param.Tstep;  % use t-P:t-1 data to predict time series vaule at time t+6
K0=param.K0; % filter order L_0 filter
K1=param.K1; % filter order L_1 filter
K2=param.K2; % filter order L_2 filter
V_En=param.V_En; %Enable node information 0 or 1
F_En=param.F_En; %Enable edge information 0 or 1
T_En=param.T_En; %Enable triangle information 0 or 1
K0=param.K0*V_En; % Mask the signal if V_En =0
K1=param.K1*F_En; % Mask the signal if F_En =0
K2=param.K2*T_En; % % Mask the signal if T_En =0
FeatureNormlzn=param.FeatureNormlzn; %enable feature normalization 0 or 1
b=param.b; %normalization parameter


L0=Hodge.L0; %L_0 Laplacian
L1_lwr=Hodge.L1_lwr; %L1_lower Laplacian
L1_upr=Hodge.L1_upr; %L1_upper Laplacian
L2=Hodge.L2; %L2 Laplacian
LL0=Hodge.LL0;  %L_0 Laplacian
LL1=Hodge.LL1;   %L_1 Laplacian
LL2=Hodge.LL2;  %L_2 Laplacian
B1=Hodge.B1; % B1 incidence matrix
B2=Hodge.B2; % B2 incidence matrix
N0=Hodge.N0; %Number of nodes
N1=Hodge.N1; %Number of edges
N2=Hodge.N2; %Number of triangles 

% Initialize bias
bias=1
if bias ==1
    bias_V=1*ones(N0,1);
    bias_F=1*ones(N1,1);
    bias_T=1*ones(N2,1);
else
    bias_V=[];
    bias_F=[];
    bias_T=[];
end
param.bias_V=bias_V;
param.bias_F=bias_F;
param.bias_T=bias_T;

%Initialization of Filter Coefficents
V_theta_old=zeros(P*K0*(K0+K1)+bias,1);
F_theta_old=zeros(P*K1*(K0+K1+K2)+bias,1);
T_theta_old=zeros(P*K2*(K1+K2)+bias,1);
%Initialization of PHI
V_PHI_old=zeros(P*K0*(K0+K1)+bias,P*K0*(K0+K1)+bias);
F_PHI_old=zeros(P*K1*(K0+K1+K2)+bias,P*K1*(K0+K1+K2)+bias);
T_PHI_old=zeros(P*K2*(K1+K2)+bias,P*K2*(K1+K2)+bias);
%Initialization r 
V_r_old=zeros(P*K0*(K0+K1)+bias,1);
F_r_old=zeros(P*K1*(K0+K1+K2)+bias,1);
T_r_old=zeros(P*K2*(K1+K2)+bias,1);


norm_scale_V=0;
norm_scale_F=0;
norm_scale_T=0;
% algorithm run in online way t=P+1 to T
for t=P+1:T
    param.t=t;
    sig_V=signal_node(:,t-P:t-1); %P time lagged values of vertex signa;
    sig_F=signal_edge(:,t-P:t-1); %P time lagged values of edge siganl
    sig_T=signal_tri(:,t-P:t-1);  %P time lagged values of triangle signal
    
    %  Estiamtion of vertex signal uses the time lagged values of vertex and edge signel
    V_00=Feature_Gen(LL0,K0,eye(N0),LL0,K0,sig_V,param);% vertex features for vertex processing
    F_01=Feature_Gen(LL0,K0,B1,LL1,K1,sig_F,param); % edge features for vertex processing
    SV=Feature_Collect(V_00,F_01,[],bias_V,param);
    if FeatureNormlzn==1 %normalzie feature of FeatureNormlzn=1
        if V_En~=0
            SV_n=sum(SV.^2);
            SV_n(SV_n==0)=0.001;
            varV=sum(V_00(:,1).^2);
            norm_scale_V=(1-b)*norm_scale_V+b*sqrt(varV);
            SV=SV./sqrt(SV_n)*norm_scale_V;
        else
            SV=0*SV;
        end
    end
    
    % %  Estiamtion of edge signal uses the time lagged values of vertex, edge  and triangle signals
    V_10=Feature_Gen(LL1,K1,B1',LL0,K0,sig_V,param);% vertex features for edge processing
    F_11=Feature_Gen(LL1,K1,eye(N1),LL1,K1,sig_F,param); % edge features for edge processing
    T_12=Feature_Gen(LL1,K1,B2,LL2,K2,sig_T,param); % triangle features for edge processing
    SF=Feature_Collect(V_10,F_11,T_12,bias_F,param);
    
    if FeatureNormlzn==1 %normalzie feature of FeatureNormlzn=1
        if F_En~=0
            SF_n=sum(SF.^2);
            SF_n(SF_n==0)=0.001;
            varF=sum(F_11(:,1).^2);
            norm_scale_F=(1-b)*norm_scale_F+b*sqrt(varF);
            SF=SF./sqrt(SF_n)*norm_scale_F;
        else
           SF=0*SF; 
        end
    end
    
     % %  Estiamtion of triangle signal uses the time lagged values of  edge  and triangle signals
    F_21=Feature_Gen(LL2,K2,B2',LL1,K1,sig_F,param); % edge features for triangle processing
    T_22=Feature_Gen(LL2,K2,eye(N2),LL2,K2,sig_T,param); % triangle features for triangle processing
    ST=Feature_Collect([],F_21,T_22,bias_T,param);
    if FeatureNormlzn==1 %normalzie feature of FeatureNormlzn=1
        if T_En~=0
            ST_n=sum(ST.^2);
            ST_n(ST_n==0)=0.001;
            varT=sum(T_22(:,1).^2);
            norm_scale_T=(1-b)*norm_scale_T+b*sqrt(varT);
            ST=ST./sqrt(ST_n)*norm_scale_T;
        else
            ST=0*ST; 
        end
    end
    
    
    %Online Estimation: Vertex
    L_v=mu_0*L0;
    [V_theta_update,V_PHI_update,V_r_update,maxeig]=SCVAR_OnlineEstimator(SV,L_v,V_PHI_old,V_r_old,V_theta_old,signal_node(:,t),param);
    V_theta_old=V_theta_update;
    V_PHI_old=V_PHI_update;
    V_r_old=V_r_update;
    V_theta_store(:,t)=V_theta_old;
    maxeig_V(t)=maxeig;
    
    %Online Estimation: Edge
    L_f=mu_1l*L1_lwr+mu_1u*L1_upr;
    [F_theta_update,F_PHI_update,F_r_update,maxeig]=SCVAR_OnlineEstimator(SF,L_f,F_PHI_old,F_r_old,F_theta_old,signal_edge(:,t),param);
    F_theta_old=F_theta_update;
    F_PHI_old=F_PHI_update;
    F_r_old=F_r_update;
    F_theta_store(:,t)=F_theta_old;
    maxeig_F(t)=maxeig;
    
    %Online Estimation: Triangle
    L_t=mu_2*L2;
    [T_theta_update,T_PHI_update,T_r_update,maxeig]=SCVAR_OnlineEstimator(ST,L_t,T_PHI_old,T_r_old,T_theta_old,signal_tri(:,t),param);
    T_theta_old=T_theta_update;
    T_PHI_old=T_PHI_update;
    T_r_old=T_r_update;
    T_theta_store(:,t)=T_theta_old;
    maxeig_T(t)=maxeig;
    
    %Intialize values for Preditction
    Logs4Predn.signal_node=signal_node;
    Logs4Predn.signal_edge=signal_edge;
    Logs4Predn.signal_tri=signal_tri;
    Logs4Predn.V_PHI_old=V_PHI_old;
    Logs4Predn.F_PHI_old=F_PHI_old;
    Logs4Predn.T_PHI_old=T_PHI_old;
    Logs4Predn.V_r_old=V_r_old;
    Logs4Predn.F_r_old=F_r_old;
    Logs4Predn.T_r_old=T_r_old;
    Logs4Predn.V_theta_old=V_theta_old;
    Logs4Predn.F_theta_old=F_theta_old;
    Logs4Predn.T_theta_old=T_theta_old;
    %  Tstep ahaed prediction of node,edge, triangle signal.
    [V_pred,F_pred,T_pred]=SCVAR_Forecast(Logs4Predn,Hodge,param);
    signal_node_pred_buffer(:,:,t+1)=V_pred;
    signal_edge_pred_buffer(:,:,t+1)=F_pred;
    signal_tri_pred_buffer(:,:,t+1)=T_pred;
    
end
% figure
% imagesc(V_theta_store)





signal_node_pred=squeeze(signal_node_pred_buffer(:,Tstep,:));    
signal_edge_pred=squeeze(signal_edge_pred_buffer(:,Tstep,:));   
signal_tri_pred=squeeze(signal_tri_pred_buffer(:,Tstep,:));   
nmse_node=mean(CompNMSE(signal_node,signal_node_pred(:,1:end-1)));
 nmse_edge=mean(CompNMSE(signal_edge,signal_edge_pred(:,1:end-1)));
nmse_tri=mean(CompNMSE(signal_tri,signal_tri_pred(:,1:end-1)));
end

% PLOT_SCVAR
% plot(nmse_node_1');
% hold on
% plot(nmse_node_2');
% legend("1","2")

