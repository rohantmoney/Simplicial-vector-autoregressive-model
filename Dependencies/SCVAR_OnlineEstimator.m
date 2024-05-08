         function [theta_update,PHI_update,r_update,maxeig]=SCVAR_OnlineEstimator(S,L,PHI_old,r_old,theta_old,x_t,param)
gamma=param.gamma;
delta=param.delta;
LassoEn=param.LassoEn;
% mu1=param.mu1;
lambda=param.lambda;
% eta=param.eta;
N=size(L,1);

PHI_update=gamma*PHI_old+delta*S'*(eye(N)+L)*S;
r_update=gamma*r_old+delta*S'*x_t;
%   eta=2*10^-7;
maxeig=max(eigs(PHI_update));
 eta=2/(0.0001+max(eigs(PHI_update)));% CHECK!!!!!
% max(eigs(PHI_update))
theta_update1=theta_old-eta*(PHI_update*theta_old-r_update+lambda*theta_old); %Gradient Descent
if LassoEn==1
   theta_update=theta_update1.*max(0.00,(1-eta*lambda./abs(theta_update1)));
else
   theta_update=theta_update1;  
end
end