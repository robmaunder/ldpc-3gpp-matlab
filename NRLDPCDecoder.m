classdef NRLDPCDecoder < matlab.System

    properties(Nontunable)
        BG = 1; % Default value
        Z = 2; % Default value
        iterations = 50; % Default value
    end
    
    properties(Nontunable, SetAccess = private)
        K
        N
        H
   end
    
    properties(Hidden, Nontunable)
        hDec
    end

    methods
           function obj = NRLDPCDecoder(varargin)
              setProperties(obj,nargin,varargin{:});
           end
    end
    
    
    methods(Access = protected)
        
        function setupImpl(obj)
            obj.H=get_pcm(get_3gpp_base_graph(obj.BG,get_3gpp_set_index(obj.Z)),obj.Z);
            obj.K = size(obj.H,2) - size(obj.H,1);            
            obj.N = size(obj.H,2);
            obj.hDec = comm.LDPCDecoder('ParityCheckMatrix',obj.H,'MaximumIterationCount',obj.iterations,'IterationTerminationCondition','Parity check satisfied');
        end

        function out = stepImpl(obj, in)
            out = step(obj.hDec, in);
        end
                  
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end
