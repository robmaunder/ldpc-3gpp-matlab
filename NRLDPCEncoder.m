classdef NRLDPCEncoder < matlab.System
    
    properties(Nontunable)
        BG = 2; % Default value
        L = 24; % Default value
        K_prime = 20; % Default value
        N_ref = 100; % Default value - swap for TBS_LBRM later
        I_LBRM = 0; % Default value
        rv_id = 0; % Default value
        E = 100; % Default value - E might be zero for some code blocks - think about this later
        Q_m = 1; % Default value
    end
    
    properties(Access = private, Hidden)
        hCRCGenerator
        hLDPCEncoder
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
        function obj = NRLDPCEncoder(varargin)
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
                obj.hCRCGenerator = comm.CRCGenerator('Polynomial',obj.CRCGeneratorPolynomial);
            end
            obj.hLDPCEncoder = comm.LDPCEncoder('ParityCheckMatrix',obj.ParityCheckMatrix);
        end
        
        function f = stepImpl(obj, b)
            c = append_CRC_and_padding(obj, b);
            d = LDPC_coding(obj, c);
            e = bit_selection(obj, d);
            f = bit_interleaving(obj, e);
        end
        
        function c = append_CRC_and_padding(obj, b)
            if length(b) ~= obj.K_prime - obj.L
                error('ldpc_3gpp_matlab:Error','Length of b should be K_prime-L.');
            end
            
            c = zeros(obj.K, 1);
            
            s = 0;
            for k = 0:obj.K_prime-obj.L-1
                c(k+1) = b(s+1);
                s = s + 1;
            end
            if obj.L == 24 % C>1
                bp = step(obj.hCRCGenerator, b);
                p = bp(obj.K_prime-obj.L+1:obj.K_prime);
                for k = obj.K_prime-obj.L:obj.K_prime-1
                    c(k+1) = p(k+obj.L-obj.K_prime+1);
                end
            end
            for k = obj.K_prime:obj.K-1
                c(k+1) = NaN;
            end
        end
        
        % Implements Section 5.3.2 of TS38.212
        function d = LDPC_coding(obj, c)
            if length(c) ~= obj.K
                error('ldpc_3gpp_matlab:Error','Length of c should be K.');
            end
            
            d = zeros(obj.N,1);
            
            % Not sure about what to do if there are NaNs within the first
            % 2*obj.Z_c elements of c. The following code (adapted from
            % TS38.212) does not set these to 0.
            for k = 2*obj.Z_c:obj.K-1
                if ~isnan(c(k+1))
                    d(k-2*obj.Z_c+1) = c(k+1);
                else
                    c(k+1) = 0;
                    d(k-2*obj.Z_c+1) = NaN;
                end
            end
            
            cw = step(obj.hLDPCEncoder, c);
            w = cw(obj.K+1:obj.N+2*obj.Z_c);
            
            for k = obj.K:obj.N+2*obj.Z_c-1
                d(k-2*obj.Z_c+1) = w(k-obj.K+1);
            end
        end
        
        % Implements Section 5.4.2.1 of TS38.212
        function e = bit_selection(obj, d)
            if length(d) ~= obj.N
                error('ldpc_3gpp_matlab:Error','Length of d should be N.');
            end
            
            e = zeros(obj.E,1);
            
            k = 0;
            j = 0;
            while k < obj.E
                if ~isnan(d(mod(obj.k_0 + j, obj.N_cb)+1))
                    e(k+1) = d(mod(obj.k_0 + j, obj.N_cb)+1);
                    k = k+1;
                end
                j = j+1;
            end
        end
        
        % Implements Section 5.4.2.2 of TS38.212
        function f = bit_interleaving(obj, e)
            if length(e) ~= obj.E
                error('ldpc_3gpp_matlab:Error','Length of e should be E.');
            end
            
            f = zeros(obj.E,1);
            
            for j = 0:obj.E/obj.Q_m-1
                for i = 0:obj.Q_m-1
                    f(i+j*obj.Q_m+1) = e(i*obj.E/obj.Q_m+j+1);
                end
            end
        end
        
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
        
    end
end
