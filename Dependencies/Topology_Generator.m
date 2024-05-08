function Hodge=Topology_Generator(param)
IncidDisplay=param.IncidDisplay;
DisplayGraph=param.DisplayGraph;
HodgeNormlzn=param.HodgeNormlzn;
edg=param.edg;
edg_index=param.edg_index;
% V_En=param.V_En;
% F_En=param.F_En;
% T_En=param.T_En;
K0=param.K0;
K1=param.K1;
K2=param.K2;
% 
% edg=[1,2;2,3;3,1;4,5;5,6;6,4];
% edg_index=[1,2,3,4,5,6];
% edg=[1,2;2,3;3,4;1,4;1,3;2,4];
% edg_index=[1,2,3,4,5,6];

[edg_sorted,edge_sort_idx]=sortrows(edg); % Lexigocraphic sorting
edg_index_sorted=edg_index(edge_sort_idx);


if DisplayGraph==1
%     figure 
    DispGraph(edg_sorted,edg_index_sorted,"Synthetic WDN")
end

[B1,B2]=Compute_Incidence(edg,param);



Hodge.B1=B1;
Hodge.B2=B2;
N0=size(B1,1); %no of nodes
N1=size(B1,2); %no of edges
N2=size(B2,2); %no of triangles
Hodge.N0=N0;
Hodge.N1=N1;
Hodge.N2=N2;

node_index_sorted=1:N0;
tringle_index_sorted=1:N2; % TAKE CARE while reding the data

% signal_node=ones(N0,T);
% signal_edge_pre=2*ones(N1,T);
% signal_tri=3*ones(N2,T);
% signal_edge=signal_edge_pre(edge_sort_idx,:);


L0=B1*B1';
L1_lwr=B1'*B1; %Lower Hodge Laplacian (Order-1)
L1_upr=B2*B2'; %Upper Hodge Laplacian (Order-1)
L2=B2'*B2;


if HodgeNormlzn==1 
    [U,S,V]=svd(L0);
    L01=U*(S/max(S(:)))*V';
    [U,S,V]=svd(L1_lwr);
    L1_lwr1=U*(S/max(S(:)))*V';
    [U,S,V]=svd(L1_upr);
    L1_upr1=U*(S/max(S(:)))*V';
    [U,S,V]=svd(L2);
    L21=U*(S/max(S(:)))*V';
end
    
Hodge.L0=L0;
Hodge.L1_lwr=L1_lwr;
Hodge.L1_upr=L1_upr;
Hodge.L2=L2;    

% *******TAKE CARE OF EDGE SORTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!

LL0=LL_Gen(L0,K0,param);
LL1_lwr=LL_Gen(L1_lwr,K1/2,param);
LL1_upr=LL_Gen(L1_upr,K1/2,param);
LL1=[LL1_lwr,LL1_upr];
LL2=LL_Gen(L2,K2,param);
Hodge.LL0=LL0;
Hodge.LL1_lwr=LL1_lwr;
Hodge.LL1_upr=LL1_upr;
Hodge.LL1=LL1;
Hodge.LL2=LL2;
end