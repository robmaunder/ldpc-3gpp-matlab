classdef NRLDPCDecoder < matlab.System
    
    properties(Nontunable)
        BG = 2; % Default value
        Z_c = 2; % Default value
        K_prime = 20; % Default value
        iterations = 50; % Default value
    end
    
    properties(Access = private, Hidden)
        hLDPCDecoder
    end
    
    properties(Dependent, SetAccess = private)
        SetIndex
        BaseGraph
        ParityCheckMatrix
        K
        N
    end
    
    methods
        function obj = NRLDPCDecoder(varargin)
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
            try
                obj.hLDPCDecoder = comm.gpu.LDPCDecoder('ParityCheckMatrix',obj.ParityCheckMatrix,'MaximumIterationCount',obj.iterations,'IterationTerminationCondition','Parity check satisfied');
            catch
                obj.hLDPCDecoder = comm.LDPCDecoder('ParityCheckMatrix',obj.ParityCheckMatrix,'MaximumIterationCount',obj.iterations,'IterationTerminationCondition','Parity check satisfied');
            end
        end
        
        function c_hat = stepImpl(obj, d_tilde)
            c_hat = LDPC_coding(obj, d_tilde);
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
            c_tilde = cw_tilde(1:obj.K);
            cw_tilde(isnan(cw_tilde)) = inf;
            c_hat = double(step(obj.hLDPCDecoder, cw_tilde));          
            c_hat(isnan(c_tilde)) = NaN;
         end
        
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
        
    end
end
