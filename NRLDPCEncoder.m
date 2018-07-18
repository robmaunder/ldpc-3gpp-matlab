classdef NRLDPCEncoder < matlab.System

    properties(Nontunable)
        BG = 1; % Default value
        Z = 2; % Default value
    end
    
    properties(Nontunable, SetAccess = private)
        K
        N
        H
   end
    
    properties(Hidden, Nontunable)
        hEnc
    end

    methods
           function obj = NRLDPCEncoder(varargin)
              setProperties(obj,nargin,varargin{:});
           end
    end
    
    
    methods(Access = protected)
        
        function setupImpl(obj)
            obj.H=get_pcm(get_3gpp_base_graph(obj.BG,get_3gpp_set_index(obj.Z)),obj.Z);
            obj.K = size(obj.H,2) - size(obj.H,1);            
            obj.N = size(obj.H,2);
            obj.hEnc = comm.LDPCEncoder('ParityCheckMatrix',obj.H);
        end

        function out = stepImpl(obj, in)
            out = step(obj.hEnc, in);
        end
                  
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end
