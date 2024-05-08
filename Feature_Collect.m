function feature_collection=Feature_Collect(V_feature,F_feature,T_feature,bias_feature,param)
  feature_collection=[];
  V_En=param.V_En;
  F_En=param.F_En;
  T_En=param.T_En;
  
  if V_En~=0
      feature_collection=[feature_collection,V_feature];
  end
  if F_En~=0
      feature_collection=[feature_collection,F_feature];
  end    
  if T_En~=0
      feature_collection=[feature_collection,T_feature];
  end
  feature_collection=[feature_collection,bias_feature];
end