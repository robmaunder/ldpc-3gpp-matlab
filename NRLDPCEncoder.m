classdef NRLDPCEncoder < matlab.System
    
    properties(Nontunable)
        BG = 2; % Default value
        Z_c = 2; % Default value
        K_prime = 20; % Default value
    end
    
    properties(Access = private, Hidden)
        hLDPCEncoder
    end
    
    properties(Dependent, SetAccess = private)
        SetIndex
        BaseGraph
        ParityCheckMatrix
        K
        N
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
        
        function set.Z_c(obj, Z_c)
            get_3gpp_set_index(Z_c); % Test if Z_c is valid
            obj.Z_c = Z_c;
        end
        
        function set.K_prime(obj, K_prime)
            if K_prime > obj.K
                error('ldpc_3gpp_matlab:UnsupportedBlockLength','K_prime should be no greater than K.');
            end
            obj.K_prime = K_prime;
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
        
    end
    
    
    methods(Access = protected)
        
        function setupImpl(obj)
            obj.hLDPCEncoder = comm.LDPCEncoder('ParityCheckMatrix',obj.ParityCheckMatrix);
        end
        
        function d = stepImpl(obj, c)
            d = LDPC_coding(obj, c);
        end
        
        % Implements Section 5.3.2 of TS38.212
        function d = LDPC_coding(obj, c)
            if length(c) ~= obj.K
                error('ldpc_3gpp_matlab:Error','Length of c should be K.');
            end
            if ~all(isnan(c(obj.K_prime+1:end)))
                error('ldpc_3gpp_matlab:Error','c should be appended with K-K_prime NaNs.');
            end
            
            d = zeros(obj.N,1);
            
            for k = 2*obj.Z_c:obj.K-1
                if ~isnan(c(k+1))
                    d(k-2*obj.Z_c+1) = c(k+1);
                else
                    c(k+1) = 0;
                    d(k-2*obj.Z_c+1) = NaN;
                end
            end
            
            cw = step(obj.hLDPCEncoder, c);
            w = cw(obj.K+1:end);
            
            for k = obj.K:obj.N+2*obj.Z_c-1
                d(k-2*obj.Z_c+1) = w(k-obj.K+1);
            end
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
        
    end
end
