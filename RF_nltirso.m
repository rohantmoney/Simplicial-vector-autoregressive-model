classdef RF_nltirso
       properties
        NoOfRF 
        lambda 
        h_stepsize       = @(RF_ts)3/eigs(RF_ts.m_Phi, 1);
        eta
        noOfNodes
        filtOrder
        vsigma
        forgettingFactor % gamma
    end
    methods
        function RF_ts = initialize(obj,sigma2_init,m_first_samples, t_alpha)
            RF_ts = RF_nltirso_state ;
            P = obj.filtOrder;
            N = obj.noOfNodes;
            D=obj.NoOfRF;
            assert(isequal(size(m_first_samples), [P, N]), ...
                'size of m_first_samples should be P x N (obj.order -by- obj.noOfNodes')
            if isscalar(sigma2_init)
                RF_ts.m_Phi = sigma2_init*eye(P*N*2*D);
            elseif ismatrix(sigma2_init)
                RF_ts.m_Phi = sigma2_init;
            else
                error 'wrong format sigma2'
            end
            RF_ts.m_R   = zeros(P*N*2*D,N); % r_n ^= ts.m_R(:,n)
            if exist('t_alpha', 'var')
                assert(isequal(size(t_alpha), [N, N, P,2*D]), ...
                    'inconsistent size of initial parameter tensor');
            else
                t_alpha = zeros(N, N, P,2*D);
                %             t_A(:,:, 1) = 0.99*eye(N);
            end
            %         ts.m_A = ts.get_m_A(t_A); % turn it into a P*N x N matrix
            RF_ts.coeff=t_alpha;
            RF_ts.m_buffer = flipud(m_first_samples);
 
            for ii=1:obj.noOfNodes
            RF_ts.v(ii,:,:)=obj.vsigma(ii)*obj.vsigma(ii)*randn(1,obj.filtOrder,D);
            end
        end
        function [RF_ts_out, alpha, m_B] = update(obj, RF_ts_in, v_y)
            alpha=RF_ts_in.coeff;
            
            D=obj.NoOfRF;
            m_X=RF_ts_in.m_buffer';
            t=obj.filtOrder+1;
            for n1=1:obj.noOfNodes
                
                for tau=1:obj.filtOrder
                    m_X_P_pre=m_X(n1,t-tau);
                    v_vec=RF_ts_in.v(n1,tau,:);
                    
                    
                    
                    Kernal(n1,tau,1:2*D)=[reshape((sin(m_X_P_pre*v_vec)),[1,D]) reshape((cos(m_X_P_pre*v_vec)),[1,D])];
                end
            end
            
            Kernal_Size_Full=size(Kernal);
            K_vec=vec(Kernal(:,:,1:2*D));
             m_Phi = obj.forgettingFactor  * RF_ts_in.m_Phi + ...
             (1-obj.forgettingFactor) * (K_vec*K_vec');
            m_R =  obj.forgettingFactor  * RF_ts_in.m_R + ...
            (1-obj.forgettingFactor) * K_vec*v_y';
         RF_ts_out=RF_ts_in.updateBuffer(v_y);
         RF_ts_out.m_Phi = m_Phi;
         RF_ts_out.m_R   = m_R;
        %RF_ts_out.eig_value_phi=max(eig((m_Phi )));
         step_re = obj.h_stepsize(RF_ts_out);   
        obj.eta=1/step_re;
            for n1=1:obj.noOfNodes
                
                grad_vec=m_Phi*(vec(alpha(n1,:,:,1:2*D)))-m_R(:,n1);
                
                grad=reshape(grad_vec,obj.noOfNodes,obj.filtOrder,2*D);
                g=obj.eta*squeeze(alpha(n1,:,:,1:2*D))-grad;
                
                alpha(n1,:,:,1:2*D)= obj.SoftThresold(g,2*D,obj.noOfNodes,obj.eta,obj.lambda,n1,obj.filtOrder);
            end
           
            RF_ts_out.coeff=alpha;
           
        end
        function [alpha_temp] =SoftThresold(obj,g,t,noOfNodes,eta,lamda,n1,order)
            alpha_temp=zeros(size(g));
            for n2=1:noOfNodes
                for m=1:order
%                              if n2==n1
%                                               alpha_temp(n2,m,:)=g(n2,m,:)/eta;
%                     % %
%                               else
%                     %
                  % alpha_temp(n2,m,:)=g(n2,m,:)/eta;
        alpha_temp(n2,m,:)=(g(n2,m,:)/eta)*max(0,(1-lamda/norm(reshape(g(n2,m,:),1,t))));
%                              end
                end
            end
        end
    end
end