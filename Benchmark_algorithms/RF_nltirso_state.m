classdef RF_nltirso_state
    properties
        m_buffer
        coeff
        v
         m_Phi    % PN2D x PN2D matrix
    m_R      % PN2D x N  matrix % collects the N PN-vectors r_n
    t    
    eig_value_phi
    end
    methods
        function obj = TirsoState(ch_message)
            % constructor. Do not use directly. Construct a Tirso and call its
            % initialize method instead.
            errormsg = ['The correct way of constructing a TirsoStateobject'...
                ' is calling the initialize method of a Tirso object'];
            if not(exist('ch_message', 'var'))
                error(errormsg)
            else
                switch ch_message
                    case 'calling from Tirso.initialize'
                        % OK
                    otherwise
                        warning (errormsg);
                end
            end
        end
        
        function ts_out = updateBuffer(obj, v_y)
            % updateBuffer: shift the buffer one position and inserts
            % a new data vector.
            % Input:  v_y    N-vector
            % Output: ts_out object encoding new state
            %
            % Note: TirsoState is not handle; in other words:
            % updateBuffer does not modify the object over which is called
            assert( iscolumn(v_y) && length(v_y)==size(obj.m_buffer,2) );
            ts_out = obj;
            ts_out.m_buffer(2:end, :) = obj.m_buffer(1:end-1, :);
            ts_out.m_buffer(1, :) = v_y;
        end
        function v_out = predictFromBuffer(obj)
            m_X=obj.m_buffer';
            alpha=obj.coeff;
            [~,t]=size(m_X);
            t=t+1;
            [noOfNodes,~,filtOrder,D]=size(alpha);
            D=D/2;
            for n1=1:noOfNodes
                
                for tau=1:filtOrder
                    m_X_P_pre=m_X(n1,t-tau);
                    v_vec=obj.v(n1,tau,:);
                    
                    
                    
                    Kernal(n1,tau,1:2*D)=[reshape((sin(m_X_P_pre*v_vec)),[1,D]) reshape((cos(m_X_P_pre*v_vec)),[1,D])];
                end
            end
            
            Kernal_Size_Full=size(Kernal);
            K_vec=vec(Kernal(:,:,1:2*D));
            for n1=1:noOfNodes
                predt(n1)=(vec(alpha(n1,:,:,1:2*D)))'*K_vec;
                m_X(n1,t)=predt(n1);
            end
            v_out=m_X(:,t);
        end
        function m_out = predictMany(obj, m_pastValues, K_toPredict)
            % predictMany predict an arbitrary (k) number of steps ahead
            % Input: m_pastValues P x N matrix containing P past value vectors
            %                                  (oldest up)
            %        K_toPredict  scalar       number k of steps to predict
            % Output: m_out       K x N matrix containing predictions for all
            %                     time horizons
            [P, N] = size(m_pastValues);
%             assert(isequal(size(obj.m_A), [P*N, N]));
            
            ts_local = obj;
            ts_local.m_buffer = flipud(m_pastValues);
            % values stored in buffer in reverse time order
            
            m_out = zeros(K_toPredict, N);
            for k = 1:K_toPredict
                m_out(k, :) = ts_local.predictFromBuffer;
                ts_local = ts_local.updateBuffer(m_out(k,:)');
            end
        end
        function m_out = predictManyFromBuffer(obj, K_toPredict)
            % predictManyFromBuffer predicts k steps ahead taking the current
            % contents of the buffer as the past data
            % Input: K_toPredict  scalar       number k of steps to predict
            
            m_out = obj.predictMany(flipud(obj.m_buffer), K_toPredict);
        end
        
    end
end