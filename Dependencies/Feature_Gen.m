function feature=Feature_Gen(L_f,K_f,B,L_s,K_s,sig,param)
  P=param.P;
  feature=[];
  if K_f*K_s~=0
      for p=1:P
          sig_p=sig(:,P+1-p);
          trans_conv= B*L_s*kron(eye(K_s),sig_p);
          feature_p=L_f*kron(eye(K_f),trans_conv);
          feature =[feature, feature_p];
      end
  end
end