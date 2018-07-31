%NRLDPCDECODER Decoder for 3GPP New Radio LDPC code
%   LDPCDEC = NRLDPCDECODER creates a 3GPP New Radio LDPC decoder system
%   object, LDPCDEC. Default values are assumed for all properties, which
%   are inherited from the NRLDPC base class.
%
%   LDPCDEC = NRLDPCDECODER(Name,Value) creates a 3GPP New Radio LDPC
%   decoder system object, LDPCDEC, with the specified property Name set to
%   the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Copyright © 2018 Robert G. Maunder. This program is free software: you 
%   can redistribute it and/or modify it under the terms of the GNU General 
%   Public License as published by the Free Software Foundation, either 
%   version 3 of the License, or (at your option) any later version. This 
%   program is distributed in the hope that it will be useful, but WITHOUT 
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
%   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
%   for more details.
classdef NRLDPCDecoder < NRLDPC
    
    properties(Nontunable)
        %I_HARQ Enable Hybrid Automatic Repeat reQuest (HARQ)
        %   Specifies whether or not successive blocks of input LLRs are assumed to
        %   be retransmissions of the same information block. When I_HARQ = 0, it
        %   is assumed that successive blocks of input LLRs correspond to different
        %   information blocks. In this case, the blocks of input LLRs are decoded
        %   independently. When I_HARQ ~= 0, the successive blocks of input LLRs
        %   are assumed to be retransmissions of the same information block. In
        %   this case the successive blocks of input LLRs are accumulated in an
        %   internal buffer, before they are processed by the LDPC decoder core.
        %   The internal buffer may be reset using the reset method, when the
        %   retransmission of a particular information block is completed, allowing
        %   the transmission of another information block.
        I_HARQ = 0; % Default value
    end
    
    properties
        %ITERATIONS Number of decoding iterations
        %   Specifies the number of flooding decoding iterations performed during
        %   LDPC decoding.
        iterations = 50; % Default value        
    end
    
    properties(Access = private, Hidden)
        %HCRCDETECTOR Cyclic Redundancy Check (CRC) detector
        %   A COMM.CRCDETECTOR used to check the CRC.  
        %
        %   See also COMM.CRCDETECTOR
        hCRCDetector

        %HLDPCDECODER Low Density Parity Check (LDPC) decoder
        %   A COMM.LDPCDECODER used to perform the LDPC decoding.  
        %
        %   See also COMM.LDPCDECODER        
        hLDPCDecoder
    end
    
    properties(DiscreteState)
        %BUFFER Internal buffer used to accumulate input LLRs during HARQ
        %   When I_HARQ ~= 0, successive blocks of input LLRs are assumed to be
        %   retransmissions of the same information block. In this case the
        %   successive blocks of input LLRs are accumulated in this internal
        %   buffer, before they are processed by the LDPC decoder core. The
        %   internal buffer may be reset using the reset method, when the
        %   retransmission of a particular information block is completed, allowing
        %   the transmission of another information block.
        buffer
    end
        
    methods
        % Constructor allowing properties to be set according to e.g.
        % a = NRLDPCDecoder('BG',1,'K_prime_minus_L',20,'E',132);
        function obj = NRLDPCDecoder(varargin)
            setProperties(obj,nargin,varargin{:});
        end        
    end
    
    methods(Access = protected)
        
        function setupImpl(obj)
            obj.hCRCDetector = comm.CRCDetector('Polynomial',obj.CRCPolynomial);
%             try
%                 obj.hLDPCDecoder = comm.gpu.LDPCDecoder('ParityCheckMatrix',obj.H,'MaximumIterationCount',obj.iterations,'IterationTerminationCondition','Parity check satisfied');
%             catch
                obj.hLDPCDecoder = comm.LDPCDecoder('ParityCheckMatrix',obj.H,'MaximumIterationCount',obj.iterations,'IterationTerminationCondition','Parity check satisfied');
%             end
            obj.buffer = zeros(obj.N,1);
            
        end
        
        function b_hat = stepImpl(obj, f_tilde)
            e_tilde = bit_interleaving(obj, f_tilde);
            d_tilde = bit_selection(obj, e_tilde);
            if obj.I_HARQ == 0
                c_hat = LDPC_coding(obj, d_tilde);
            else
                obj.buffer = obj.buffer + d_tilde;            
                c_hat = LDPC_coding(obj, obj.buffer);
            end
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
                for k = obj.K_prime_minus_L:obj.K_prime-1
                    p_hat(k+obj.L-obj.K_prime+1) = c_hat(k+1);
                end
                bp_hat = [b_hat; p_hat];                
                [~,err] = step(obj.hCRCDetector,bp_hat);
                if err
                    b_hat = [];
                end
            
        end
                
        function resetImpl(obj)
            obj.buffer = zeros(obj.N,1);
        end
        
    end
end
