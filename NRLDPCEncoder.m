classdef NRLDPCEncoder < matlab.System
    
    properties(Nontunable)
        BG = 1; % Default value
        Z = 2; % Default value
    end
    
    properties(SetAccess = private)
        hEnc
    end
    
    properties(Dependent)
        K
        N
    end
    
    methods
        function obj = NRLDPCEncoder(varargin)
            setProperties(obj,nargin,varargin{:});
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
            H=get_pcm(get_3gpp_base_graph(obj.BG,get_3gpp_set_index(obj.Z)),obj.Z);
            obj.hEnc = comm.LDPCEncoder('ParityCheckMatrix',H);
        end
        
        function out = stepImpl(obj, in)
            out = step(obj.hEnc, in);
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end
