classdef NRLDPCDecoder < matlab.System
    
    properties(Nontunable)
        BG = 1; % Default value
        Z = 2; % Default value
        iterations = 50; % Default value
    end
 
    properties(Access = private, Hidden)
         hDec
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
            if BG ~= 1 && BG ~= 2
                error('ldpc_3gpp_matlab:UnsupportedBaseGraph','Valid values of BG are 1 and 2.');
            end
            obj.BG = BG;
        end
        
        function set.Z(obj, Z)
            get_3gpp_set_index(Z); % Test if Z is valid
            obj.Z = Z;
        end
        
        function SetIndex = get.SetIndex(obj)
            SetIndex = get_3gpp_set_index(obj.Z);
        end
        
        function BaseGraph = get.BaseGraph(obj)
            BaseGraph = get_3gpp_base_graph(obj.BG,obj.SetIndex);
        end
        
        function ParityCheckMatrix = get.ParityCheckMatrix(obj)
            ParityCheckMatrix = get_pcm(obj.BaseGraph,obj.Z);
        end        
        
        function K = get.K(obj)
            if obj.BG == 1
                K = obj.Z*22;
            elseif obj.BG == 2
                K = obj.Z*10;
            else
                error('ldpc_3gpp_matlab:UnsupportedBaseGraph','Valid values of BG are 1 and 2.');
            end
        end
        
        function N = get.N(obj)
            if obj.BG == 1
                N = obj.Z*68;
            elseif obj.BG == 2
                N = obj.Z*52;
            else
                error('ldpc_3gpp_matlab:UnsupportedBaseGraph','Valid values of BG are 1 and 2.');
            end
        end       
    end
    
    methods(Access = protected)
        
        function setupImpl(obj)
            try
                obj.hDec = comm.gpu.LDPCDecoder('ParityCheckMatrix',obj.ParityCheckMatrix,'MaximumIterationCount',obj.iterations,'IterationTerminationCondition','Parity check satisfied');
            catch
                obj.hDec = comm.LDPCDecoder('ParityCheckMatrix',obj.ParityCheckMatrix,'MaximumIterationCount',obj.iterations,'IterationTerminationCondition','Parity check satisfied');
            end
        end
        
        function out = stepImpl(obj, in)
            out = step(obj.hDec, in);
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
        
    end
end
