classdef NRLDPC < matlab.System
    
    properties(Nontunable)
        BG = 1; % Default value
        CRC = 'CRC24B'; % Default value
        K_prime_minus_L = 20; % Default value
        N_ref = 132; % Default value - swap for TBS_LBRM later
        I_LBRM = 0; % Default value
    end
    
    properties
        rv_id = 0; % Default value
        E = 132; % Default value - E might be zero for some code blocks - think about this later
        Q_m = 1; % Default value
    end
        
    properties(Dependent, SetAccess = protected)
        CRCPolynomial
        L
        K_prime
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
        function obj = NRLDPC(varargin)
            setProperties(obj,nargin,varargin{:});
        end
        
        function set.BG(obj, BG)
            if BG < 1 || BG > 2
                error('ldpc_3gpp_matlab:UnsupportedBaseGraph','Valid values of BG are 1 and 2.');
            end
            obj.BG = BG;
        end
        
        function set.CRC(obj, CRC)
            get_3gpp_crc_polynomial(CRC);
            obj.CRC = CRC;
        end
        
        function set.K_prime_minus_L(obj, K_prime_minus_L)
            if K_prime_minus_L < 0
                error('ldpc_3gpp_matlab:UnsupportedBlockLength','K_prime_minus_L should not be negative.');
            end
            obj.K_prime_minus_L = K_prime_minus_L;
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
        
        function CRCPolynomial = get.CRCPolynomial(obj)
            CRCPolynomial = get_3gpp_crc_polynomial(obj.CRC);
        end

        function L = get.L(obj)
            [~,L] = get_3gpp_crc_polynomial(obj.CRC);
        end

        function K_prime = get.K_prime(obj)
            K_prime = obj.K_prime_minus_L + obj.L;
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

        end
        
        function stepImpl(obj)

        end
       
        function resetImpl(obj)

        end
        
    end
end
