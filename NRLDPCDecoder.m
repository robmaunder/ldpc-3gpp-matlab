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

        %ITERATIONS Number of decoding iterations
        %   Specifies the number of flooding decoding iterations performed during
        %   LDPC decoding.
        MaximumIterationCount = 50; % Default value    
                
        DecisionMethod = 'Hard decision';
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
%                 obj.hLDPCDecoder = comm.gpu.LDPCDecoder('ParityCheckMatrix',obj.H,'MaximumIterationCount',obj.MaximumIterationCount,'IterationTerminationCondition','Parity check satisfied', 'DecisionMethod', obj.DecisionMethod, 'OutputValue', 'Whole codeword');
%             catch
                obj.hLDPCDecoder = comm.LDPCDecoder('ParityCheckMatrix',obj.H,'MaximumIterationCount',obj.MaximumIterationCount,'IterationTerminationCondition','Parity check satisfied', 'DecisionMethod', obj.DecisionMethod, 'OutputValue', 'Whole codeword');
%             end
            obj.buffer = zeros(obj.N,1);
            
        end
        
        function b_hat = stepImpl(obj, f_tilde)
            e_tilde = bit_deinterleaving(obj, f_tilde);
            d_tilde = bit_deselection(obj, e_tilde);
            if obj.I_HARQ == 0
                c_hat = LDPC_decoding(obj, d_tilde);
            else
                obj.buffer = obj.buffer + d_tilde;            
                c_hat = LDPC_decoding(obj, obj.buffer);
            end
            b_hat = remove_CRC_and_padding(obj, c_hat);
        end
        
        % Implements Section 5.4.2.2 of TS38.212
        function e_tilde = bit_deinterleaving(obj, f_tilde)
            if size(f_tilde,1) ~= obj.E || size(f_tilde,2) ~= 1
                error('ldpc_3gpp_matlab:Error','f_tilde should be a column vector of length E.');
            end
            
            e_tilde = zeros(obj.E,1);
            
            for j = 0:obj.E/obj.Q_m-1
                for i = 0:obj.Q_m-1
                    e_tilde(i*obj.E/obj.Q_m+j+1) = f_tilde(i+j*obj.Q_m+1);
                end
            end
        end
        
        % Implements Section 5.4.2.1 of TS38.212
        function d_tilde = bit_deselection(obj, e_tilde)
            if size(e_tilde,1) ~= obj.E || size(e_tilde,2) ~= 1
                error('ldpc_3gpp_matlab:Error','e_tilde should be a column vector of length E.');
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
        function [c_tilde_p, d_tilde_p] = LDPC_decoding(obj, d_tilde_a, c_tilde_a)
            if size(d_tilde_a,1) ~= obj.N || size(d_tilde_a,2) ~= 1
                error('ldpc_3gpp_matlab:Error','d_tilde_a should be a column vector of length N.');
            end
            if nargin == 3
                if ~strcmp(obj.DecisionMethod,'Soft decision')
                    error('ldpc_3gpp_matlab:Error','c_tilde_a input is only supported for DecisionMethod of ''Soft decision''');
                end
                if size(c_tilde_a,1) ~= obj.K || size(c_tilde_a,2) ~= 1
                    error('ldpc_3gpp_matlab:Error','c_tilde_a should be a column vector of length K.');
                end
                if ~isequal(isnan(c_tilde_a(2*obj.Z_c+1:end)),isnan(d_tilde_a(1:obj.K - 2*obj.Z_c)))
                    error('ldpc_3gpp_matlab:Error','c_tilde_a and d_tilde_a should have padding bits in corresponding positions.');
                end
            end
                            
            cw_tilde_a = [zeros(2*obj.Z_c,1); d_tilde_a];
            padding_mask = isnan(cw_tilde_a);            
            if nargin == 3
                cw_tilde_a(1:obj.K) = cw_tilde_a(1:obj.K) + c_tilde_a;
            end            
            cw_tilde_a(padding_mask) = inf;
            cw_tilde_p = double(step(obj.hLDPCDecoder, cw_tilde_a));
            cw_tilde_p(padding_mask) = NaN;
            c_tilde_p = cw_tilde_p(1:obj.K);
            d_tilde_p = cw_tilde_p(2*obj.Z_c+1:end);
        end
        
        function b_hat = remove_CRC_and_padding(obj, c_hat)
            if size(c_hat,1) ~= obj.K || size(c_hat,2) ~= 1
                error('ldpc_3gpp_matlab:Error','c_hat should be a column vector of length K.');
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
