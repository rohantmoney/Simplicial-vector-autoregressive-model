function LL=LL_Gen(L,K,param)
HodgeNormlzn=param.HodgeNormlzn;
LL=[];
if K~=0
    for i=0:K-1
        LL=[LL,L^i];
    end
    if HodgeNormlzn==1
        [U,S,V]=svd(LL);
        LL=U*S/max(S(:))*V';
    end
end
    
end