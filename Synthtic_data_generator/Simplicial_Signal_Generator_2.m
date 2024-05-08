function [signal_node,signal_edge,signal_tri]=Simplicial_Signal_Generator(Hodge,param)
% L0=Hodge.L0;
L1_lwr=Hodge.L1_lwr;
% L1_upr=Hodge.L1_upr;
% L2=Hodge.L2;
B1=Hodge.B1;
B2=Hodge.B2;
N0=Hodge.N0;
N1=Hodge.N1;
N2=Hodge.N2;
T=param.T;
% P=1;

N=N0+N1+N2;
% sigma=0.1;
CovMx=0.5*ones(N,N);
CovMx=CovMx-diag(diag(CovMx))+diag(ones(N,1));
mu=zeros(N,1);
noise_sig=mvnrnd(mu,CovMx,T)';
% noise_sig =sigma*randn(N,T);
const=0*randi(5,1);



[VAR_comb_mx1] =  VAR_MX_GEN(N);
[VAR_comb_mx2] =  VAR_MX_GEN(N);
% sig_comb(:,1)= 1*randn(N,1);
sig_comb(:,1)=1*mvnrnd(randn(N,1),CovMx,1)';
sig_comb(:,2)=1*mvnrnd(randn(N,1),CovMx,1)';
f1=1;
f2=0;
for t=3:T
    if mod(t,4000)==0
        [VAR_comb_mx1] =  VAR_MX_GEN(N);
        [VAR_comb_mx2] =  VAR_MX_GEN(N);
        swap=f1;
        f1=f2;
        f2=swap;
    end
     sig_comb(:,t)=f1*VAR_comb_mx1*sig_comb(:,t-1)+f2*VAR_comb_mx2*sig_comb(:,t-2)+noise_sig(:,t)+const; 
   % sig_comb(:,t)=VAR_comb_mx1*sig_comb(:,t-1)+noise_sig(:,t)+const;  
%    Sn=sig_comb(1:N0,t);
% Se=sig_comb(N0+1:N0+N1,t);
% St=sig_comb(N0+N1+1:end,t);


% signal_node(:,t)=1*nl_fn(Sn)+1*(B1*Se);
% signal_edge(:,t)=1*nl_fn(Se)+1*(B1'*Sn)+1*(B2*St);
% signal_tri(:,t)=1*nl_fn(St)+1*(B2'*Se);

end
Sn=sig_comb(1:N0,:);
Se=sig_comb(N0+1:N0+N1,:);
St=sig_comb(N0+N1+1:end,:);


signal_node=12*nl_fn(Sn)+1*(B1*Se);
signal_edge=12*nl_fn(Se)+1*(B1'*Sn)+1*(B2*St);
signal_tri=12*nl_fn(St)+1*(B2'*Se);


% figure
% subplot(2,1,1)
% plot(sig_comb(1:5,:)');title("Signal")
% subplot(2,1,2)
% plot(noise_sig(1:5,:)');title("Noise")
%  figure
% imagesc(VAR_comb_mx1)
grid on
function [VAR_comb_mx]=  VAR_MX_GEN(N)
    A=randn(N);
    [U,D]=eig(A);
    D2=D/max(abs(D(:)))*1.0;
    VAR_comb_mx=real(U*D2*U^-1);
%     imagesc(VAR_comb_mx);
%     grid on
    
%     D= unifrnd(0.85,0.9,N,1);
%     if mod(N,2)==0
%        sign=repmat([1,-1],1,N/2)';
%     else
%        sign=repmat([1,-1],1,(N-1)/2)'; 
%        sign=[sign;1];
%     end
%     D=D.*sign;
%     A_temp=magic(size(D,1));
%    [Q,~]=qr(A_temp);
%    VAR_comb_mx=Q'*diag(D)*Q;
end

function [VAR_comb_mx]=  change_MX(A)
    mask=randi([0,1],size(A));
    VAR_comb_mx=A.*mask;
    [U,D]=eig(VAR_comb_mx);
    D2=D/max(abs(D(:)))*1.0;
    VAR_comb_mx=real(U*D2*U^-1);
end


function y = nl_fn(x)
y =(sin(x));
end


end






% figure
% z=eig(A);
% x=0:2*pi/1000:2*pi;
% plot(real(z),imag(z),"o")
% hold on
% plot(sin(x),cos(x))
% axis equal
% grid on
% xlabel("Re(z)")
% ylabel("Im(z)")
