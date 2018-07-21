classdef NRLDPCDecoder2 < matlab.System
    
    properties(Nontunable)
        BG = 2; % Default value
        L = 24; % Default value
        K_prime = 20; % Default value
        N_ref = 100; % Default value - swap for TBS_LBRM later
        I_LBRM = 0; % Default value
        rv_id = 0; % Default value
        E = 100; % Default value - E might be zero for some code blocks - think about this later
        Q_m = 1; % Default value
        iterations = 100; % Default value
    end
    
    properties(Access = private, Hidden)
        hCRCDetector
        hLDPCDecoder
    end
    
    properties(Dependent, SetAccess = private)
        CRCGeneratorPolynomial        
        K_b
        SetIndex
        Z_c
        BaseGraph
        ParityCheckMatrix
        K
        N
        N_cb
        k_0
    end
    
    methods
        function obj = NRLDPCDecoder2(varargin)
            setProperties(obj,nargin,varargin{:});
        end
        
        function set.BG(obj, BG)
            if BG < 1 || BG > 2
                error('ldpc_3gpp_matlab:UnsupportedBaseGraph','Valid values of BG are 1 and 2.');
            end
            obj.BG = BG;
        end
        
        function set.L(obj, L)
            if L ~= 0 && L ~= 24
                error('ldpc_3gpp_matlab:UnsupportedCRCLength','Valid values of L are 0 and 24.');
            end
            obj.L = L;
        end
        
        function set.K_prime(obj, K_prime)
            if K_prime < 0
                error('ldpc_3gpp_matlab:UnsupportedBlockLength','K_prime should not be negative.');
            end
            obj.K_prime = K_prime;
        end
        
        function set.N_ref(obj, N_ref)
            if N_ref < 0
                error('ldpc_3gpp_matlab:UnsupportedBlockLength','N_ref should not be negative.');
            end
            obj.N_ref = N_ref;
        end
               
        function set.I_LBRM(obj, I_LBRM)
            obj.I_LBRM = I_LBRM;
        end
        
        function set.rv_id(obj, rv_id)
            if rv_id < 0 || rv_id > 3
                error('ldpc_3gpp_matlab:UnsupportedParameter','Valid values of rv_id are 0, 1, 2 and 3.');
            end
            obj.rv_id = rv_id;
        end
        
        function set.E(obj, E)
            if E < 0
                error('ldpc_3gpp_matlab:UnsupportedBlockLength','E should not be negative.');
            end
            obj.E = E;
        end
               
        function set.Q_m(obj, Q_m)
            if Q_m <= 0
                error('ldpc_3gpp_matlab:UnsupportedBlockLength','Q_m should be positive.');
            end
            obj.Q_m = Q_m;
        end

        function CRCGeneratorPolynomial = get.CRCGeneratorPolynomial(obj)
            if obj.L == 24
                CRCGeneratorPolynomial = 'z^24 + z^23 + z^6 + z^5 + z + 1';
            elseif obj.L == 0
                CRCGeneratorPolynomial = '';
            else
                error('ldpc_3gpp_matlab:UnsupportedCRCLength','Valid values of L are 0 and 24.');
            end            
        end
                
        function K_b = get.K_b(obj)
            if obj.BG == 1
                K_b = 22;
            elseif obj.BG == 2
                 % TS 38.212 uses B rather than K_prime for the comparisons
                 % below, but both ways give the same answer in all cases
                if obj.K_prime > 640
                    K_b = 10;
                elseif obj.K_prime > 560
                    K_b = 9;
                elseif obj.K_prime > 192
                    K_b = 8;
                else
                    K_b = 6;
                end
            else
                error('ldpc_3gpp_matlab:UnsupportedBaseGraph','BG must be 1 or 2');
            end
        end
        
        function Z_c = get.Z_c(obj)
            Z_c = get_3gpp_lifting_size(obj.K_b, obj.K_prime);
        end
               
        function SetIndex = get.SetIndex(obj)
            SetIndex = get_3gpp_set_index(obj.Z_c);
        end
        
        function BaseGraph = get.BaseGraph(obj)
            BaseGraph = get_3gpp_base_graph(obj.BG,obj.SetIndex);
        end
        
        function ParityCheckMatrix = get.ParityCheckMatrix(obj)
            ParityCheckMatrix = get_pcm(obj.BaseGraph,obj.Z_c);
            ParityCheckMatrix = [ParityCheckMatrix(:,1:obj.K_prime),ParityCheckMatrix(:,obj.K+1:end)];
        end
        
        function K = get.K(obj)
            if obj.BG == 1
                K = obj.Z_c*22;
            elseif obj.BG == 2
                K = obj.Z_c*10;
            else
                error('ldpc_3gpp_matlab:UnsupportedBaseGraph','Valid values of BG are 1 and 2.');
            end
        end
        
        function N = get.N(obj)
            if obj.BG == 1
                N = obj.Z_c*66;
            elseif obj.BG == 2
                N = obj.Z_c*50;
            else
                error('ldpc_3gpp_matlab:UnsupportedBaseGraph','Valid values of BG are 1 and 2.');
            end
        end
        
        function N_cb = get.N_cb(obj)
            if obj.I_LBRM == 0
                N_cb = obj.N;
            else
                N_cb = min(obj.N, obj.N_ref);
            end
        end
        
        function k_0 = get.k_0(obj)
            if obj.BG == 1
                if obj.rv_id == 0
                    k_0 = 0;                    
                elseif obj.rv_id == 1
                    k_0 = floor((17*obj.N_cb)/(66*obj.Z_c))*obj.Z_c;
                elseif obj.rv_id == 2
                    k_0 = floor((33*obj.N_cb)/(66*obj.Z_c))*obj.Z_c;
                elseif obj.rv_id == 3
                    k_0 = floor((56*obj.N_cb)/(66*obj.Z_c))*obj.Z_c;
                else
                    error('ldpc_3gpp_matlab:UnsupportedRedundancyVersion','Valid values of rv_id are 0, 1, 2 and 3.')
                end
            elseif obj.BG == 2
                if obj.rv_id == 0
                    k_0 = 0;                    
                elseif obj.rv_id == 1
                    k_0 = floor((13*obj.N_cb)/(50*obj.Z_c))*obj.Z_c;
                elseif obj.rv_id == 2
                    k_0 = floor((25*obj.N_cb)/(50*obj.Z_c))*obj.Z_c;
                elseif obj.rv_id == 3
                    k_0 = floor((43*obj.N_cb)/(50*obj.Z_c))*obj.Z_c;
                else
                    error('ldpc_3gpp_matlab:UnsupportedRedundancyVersion','Valid values of rv_id are 0, 1, 2 and 3.')
                end
            else
                error('ldpc_3gpp_matlab:UnsupportedBaseGraph','Valid values of BG are 1 and 2.');
            end
        end            
        
    end
    
    methods(Access = protected)
        
        function setupImpl(obj)
            if obj.L == 24
                obj.hCRCDetector = comm.CRCDetector('Polynomial',obj.CRCGeneratorPolynomial);
            end            
%             try
%                 obj.hLDPCDecoder = comm.gpu.LDPCDecoder('ParityCheckMatrix',obj.ParityCheckMatrix,'MaximumIterationCount',obj.iterations,'IterationTerminationCondition','Parity check satisfied');
%             catch
                obj.hLDPCDecoder = comm.LDPCDecoder('ParityCheckMatrix',obj.ParityCheckMatrix,'MaximumIterationCount',obj.iterations,'IterationTerminationCondition','Parity check satisfied');
%             end
        end
        
        function b_hat = stepImpl(obj, f_tilde)
            e_tilde = bit_interleaving(obj, f_tilde);
            d_tilde = bit_selection(obj, e_tilde);
            c_hat = LDPC_coding(obj, d_tilde);
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
            if ~all(isnan(d_tilde(max(obj.K_prime-2*obj.Z_c+1,1):obj.K-2*obj.Z_c)))
                error('ldpc_3gpp_matlab:Error','d_tilde should have NaNs in the appropriate places.');
            end
            
            cw_tilde = [zeros(2*obj.Z_c,1); d_tilde];
            c_hat = cw_tilde(1:obj.K);
            cw_tilde = cw_tilde(~isnan(cw_tilde));
            c_hat(~isnan(c_hat)) = double(step(obj.hLDPCDecoder, cw_tilde));          
         end
        
        function b_hat = append_CRC_and_padding(obj, c_hat)
            if length(c_hat) ~= obj.K
                error('ldpc_3gpp_matlab:Error','Length of c_hat should be K.');
            end
            
            b_hat = zeros(obj.K_prime - obj.L, 1);
            p_hat = zeros(obj.L,1);
            
            s = 0;
            for k = 0:obj.K_prime-obj.L-1
                b_hat(s+1) = c_hat(k+1);
                s = s + 1;
            end
            if obj.L == 24 % C>1
                for k = obj.K_prime-obj.L:obj.K_prime-1
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
            % Initialize / reset discrete-state properties
        end
        
    end
end
