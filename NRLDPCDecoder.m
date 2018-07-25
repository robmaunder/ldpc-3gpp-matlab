classdef NRLDPCDecoder < NRLDPC
    
    properties(Nontunable)
        I_HARQ = 0;
    end
    
    properties
        iterations = 50; % Default value        
    end
    
    properties(Access = private, Hidden)
        hCRCDetector
        hLDPCDecoder
    end
    
    properties(DiscreteState)
        d_tilde2
    end
        
    methods
        function obj = NRLDPCDecoder(varargin)
            setProperties(obj,nargin,varargin{:});
        end        
    end
    
    methods(Access = protected)
        
        function setupImpl(obj)
            if obj.L == 24
                obj.hCRCDetector = comm.CRCDetector('Polynomial',obj.CRCPolynomial);
            end
%             try
%                 obj.hLDPCDecoder = comm.gpu.LDPCDecoder('ParityCheckMatrix',obj.ParityCheckMatrix,'MaximumIterationCount',obj.iterations,'IterationTerminationCondition','Parity check satisfied');
%             catch
                obj.hLDPCDecoder = comm.LDPCDecoder('ParityCheckMatrix',obj.ParityCheckMatrix,'MaximumIterationCount',obj.iterations,'IterationTerminationCondition','Parity check satisfied');
%             end
            obj.d_tilde2 = zeros(obj.N,1);
            
        end
        
        function b_hat = stepImpl(obj, f_tilde)
            e_tilde = bit_interleaving(obj, f_tilde);
            d_tilde = bit_selection(obj, e_tilde);
            if obj.I_HARQ == 0
                obj.d_tilde2 = d_tilde;
            else
                obj.d_tilde2 = obj.d_tilde2 + d_tilde;            
            end
            c_hat = LDPC_coding(obj, obj.d_tilde2);
            b_hat = append_CRC_and_padding(obj, c_hat);
        end
        
        % Implements Section 5.4.2.2 of TS38.212
        function e_tilde = bit_interleaving(obj, f_tilde)
            if length(f_tilde) ~= obj.E
                error('ldpc_3gpp_matlab:Error','Length of f_tilde should be E.');
            end
            
            e_tilde = zeros(obj.E,1);
            
            for j = 0:obj.E/obj.Q_m-1
                for i = 0:obj.Q_m-1
                    e_tilde(i*obj.E/obj.Q_m+j+1) = f_tilde(i+j*obj.Q_m+1);
                end
            end
        end
        
        % Implements Section 5.4.2.1 of TS38.212
        function d_tilde = bit_selection(obj, e_tilde)
            if length(e_tilde) ~= obj.E
                error('ldpc_3gpp_matlab:Error','Length of e_tilde should be E.');
            end
            
            d_tilde = zeros(obj.N,1);
            d_tilde(max(obj.K_prime-2*obj.Z_c+1,1):obj.K-2*obj.Z_c) = NaN;
            
            k = 0;
            j = 0;
            while k < obj.E
                if ~isnan(d_tilde(mod(obj.k_0 + j, obj.N_cb)+1))
                    d_tilde(mod(obj.k_0 + j, obj.N_cb)+1) = d_tilde(mod(obj.k_0 + j, obj.N_cb)+1) + e_tilde(k+1);
                    k = k+1;
                end
                j = j+1;
            end
        end
        
        % Implements Section 5.3.2 of TS38.212
        function c_hat = LDPC_coding(obj, d_tilde)
            if length(d_tilde) ~= obj.N
                error('ldpc_3gpp_matlab:Error','Length of d_tilde should be N.');
            end
            
            cw_tilde = [zeros(2*obj.Z_c,1); d_tilde];
            c_tilde = cw_tilde(1:obj.K);
            cw_tilde(isnan(cw_tilde)) = inf;
            c_hat = double(step(obj.hLDPCDecoder, cw_tilde));
            c_hat(isnan(c_tilde)) = NaN;
        end
        
        function b_hat = append_CRC_and_padding(obj, c_hat)
            if length(c_hat) ~= obj.K
                error('ldpc_3gpp_matlab:Error','Length of c_hat should be K.');
            end
            
            b_hat = zeros(obj.K_prime_minus_L, 1);
            p_hat = zeros(obj.L,1);
            
            s = 0;
            for k = 0:obj.K_prime_minus_L-1
                b_hat(s+1) = c_hat(k+1);
                s = s + 1;
            end
            if obj.L == 24 % C>1
                for k = obj.K_prime_minus_L:obj.K_prime-1
                    p_hat(k+obj.L-obj.K_prime+1) = c_hat(k+1);
                end
                bp_hat = [b_hat; p_hat];                
                [~,err] = step(obj.hCRCDetector,bp_hat);
                if err
                    b_hat = [];
                end
            end
            
        end
        
        
        
        function resetImpl(obj)
            obj.d_tilde2 = zeros(obj.N,1);
        end
        
    end
end
