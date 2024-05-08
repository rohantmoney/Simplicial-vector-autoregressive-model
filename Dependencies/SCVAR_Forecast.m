
function [V_pred,F_pred,T_pred]=SCVAR_Forecast(Logs4Predn,Hodge,param)
% function [V_pred,F_pred,T_pred]=SCVAR_Forecast(signal_node,signal_edge,signal_tri,V_theta_update,F_theta_update,T_theta_update,Hodge,param)
t=param.t;
P=param.P;
V_En=param.V_En;
F_En=param.F_En;
T_En=param.T_En;
K0=param.K0*V_En;
K1=param.K1*F_En;
K2=param.K2*T_En;
mu_0=param.mu_0;
mu_1l=param.mu_1l;
mu_1u=param.mu_1u;
mu_2=param.mu_2;
lambda=param.lambda;
gamma=param.gamma;
LL0=Hodge.LL0;
LL1=Hodge.LL1;
LL2=Hodge.LL2;
B1=Hodge.B1;
B2=Hodge.B2;
L0=Hodge.L0;
L1_lwr=Hodge.L1_lwr;
L1_upr=Hodge.L1_upr;
L2=Hodge.L2;
T=param.T;
Tstep=param.Tstep;
N0=Hodge.N0;
N1=Hodge.N1;
N2=Hodge.N2;
bias_V=param.bias_V;
bias_F=param.bias_F;
bias_T=param.bias_T;
V_En=param.V_En;
F_En=param.F_En;
T_En=param.T_En;
FeatureNormlzn=param.FeatureNormlzn;
b=param.b;

signal_node=Logs4Predn.signal_node;
signal_edge=Logs4Predn.signal_edge;
signal_tri=Logs4Predn.signal_tri;
V_PHI_old=Logs4Predn.V_PHI_old;
F_PHI_old=Logs4Predn.F_PHI_old;
T_PHI_old=Logs4Predn.T_PHI_old;
V_r_old=Logs4Predn.V_r_old;
F_r_old=Logs4Predn.F_r_old;
T_r_old=Logs4Predn.T_r_old;
V_theta_old=Logs4Predn.V_theta_old;
F_theta_old=Logs4Predn.F_theta_old;
T_theta_old=Logs4Predn.T_theta_old;






    sig_V=signal_node(:,t-P+1:t);
    sig_F=signal_edge(:,t-P+1:t);
    sig_T=signal_tri(:,t-P+1:t);
    V_sig_hat=[];
    F_sig_hat=[];
    T_sig_hat=[];
    norm_scale_V=0;
    norm_scale_F=0;
    norm_scale_T=0;

for k=1:Tstep
    st_index=min(k,2);
    sig_V=[sig_V(:,st_index:end),V_sig_hat];
    sig_F=[sig_F(:,st_index:end),F_sig_hat];
    sig_T=[sig_T(:,st_index:end),T_sig_hat];

    % vertex processing features
    V_00=Feature_Gen(LL0,K0,eye(N0),LL0,K0,sig_V,param);% vertex features for vertex processing
    F_01=Feature_Gen(LL0,K0,B1,LL1,K1,sig_F,param); % edge features for vertex processing
    SV=Feature_Collect(V_00,F_01,[],bias_V,param);
    if FeatureNormlzn==1
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
 
    % edge processing features
    V_10=Feature_Gen(LL1,K1,B1',LL0,K0,sig_V,param);% vertex features for edge processing
    F_11=Feature_Gen(LL1,K1,eye(N1),LL1,K1,sig_F,param); % edge features for edge processing
    T_12=Feature_Gen(LL1,K1,B2,LL2,K2,sig_T,param); % triangle features for edge processing
    SF=Feature_Collect(V_10,F_11,T_12,bias_F,param);
    if FeatureNormlzn==1
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
    
     % triangle processing features
    F_21=Feature_Gen(LL2,K2,B2',LL1,K1,sig_F,param); % edge features for triangle processing
    T_22=Feature_Gen(LL2,K2,eye(N2),LL2,K2,sig_T,param); % triangle features for triangle processing
    ST=Feature_Collect([],F_21,T_22,bias_T,param);
    if FeatureNormlzn==1
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
    
    V_sig_hat=SV*V_theta_old;
    F_sig_hat=SF*F_theta_old;
    T_sig_hat=ST*T_theta_old;
    V_pred(:,k)=V_sig_hat;
    F_pred(:,k)=F_sig_hat;
    T_pred(:,k)=T_sig_hat;
    
    
    %Online Estimation: Vertex
    L_v=mu_0*L0;
    [V_theta_update,V_PHI_update,V_r_update,maxeig]=SCVAR_OnlineEstimator(SV,L_v,V_PHI_old,V_r_old,V_theta_old,V_sig_hat,param);
    V_theta_old=V_theta_update;
    V_PHI_old=V_PHI_update;
    V_r_old=V_r_update;
    
    %Online Estimation: Edge
    L_f=mu_1l*L1_lwr+mu_1u*L1_upr;
    [F_theta_update,F_PHI_update,F_r_update,maxeig]=SCVAR_OnlineEstimator(SF,L_f,F_PHI_old,F_r_old,F_theta_old,F_sig_hat,param);
    F_theta_old=F_theta_update;
    F_PHI_old=F_PHI_update;
    F_r_old=F_r_update;
    
    %Online Estimation: Triangle
    L_t=mu_2*L2;
    [T_theta_update,T_PHI_update,T_r_update,maxeig]=SCVAR_OnlineEstimator(ST,L_t,T_PHI_old,T_r_old,T_theta_old,T_sig_hat,param);
    T_theta_old=T_theta_update;
    T_PHI_old=T_PHI_update;
    T_r_old=T_r_update;
    

    
end




