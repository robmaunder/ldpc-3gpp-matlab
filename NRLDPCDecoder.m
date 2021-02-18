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
        %   internal d_tilde_buffer, before they are processed by the LDPC decoder core.
        %   The internal d_tilde_buffer may be reset using the reset method, when the
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
        %HTBCRCDETECTOR Cyclic Redundancy Check (CRC) detector
        %   A COMM.CRCDETECTOR used to check the transport block CRC.
        %
        %   See also COMM.CRCDETECTOR
        hTBCRCDetector
        
        %HCBCRCDETECTOR Cyclic Redundancy Check (CRC) detector
        %   A COMM.CRCDETECTOR used to check the code block CRC.
        %
        %   See also COMM.CRCDETECTOR
        hCBCRCDetector
        
        %HLDPCDECODER Low Density Parity Check (LDPC) decoder
        %   A COMM.LDPCDECODER used to perform the LDPC decoding.
        %
        %   See also COMM.LDPCDECODER
        hLDPCDecoder
    end
    
    properties(DiscreteState)
        %D_TILDE_BUFFER Internal buffer used to accumulate input LLRs during HARQ
        %   When I_HARQ ~= 0, successive blocks of input LLRs are assumed to be
        %   retransmissions of the same information block. In this case, the
        %   successive blocks of input LLRs are accumulated in this internal
        %   buffer, before they are processed by the LDPC decoder core. The
        %   internal buffer may be reset using the reset method, when the
        %   retransmission of a particular information block is completed, allowing
        %   the transmission of another information block.
        d_tilde_buffer
        
        %B_HAT_BUFFER Internal buffer used to collect successfully decoded bits during HARQ
        %   When I_HARQ ~= 0, successive blocks of input LLRs are assumed
        %   to be retransmissions of the same information block. However,
        %   some retransmissions may not contain all code blocks. In this
        %   case, the successfully decoded code blocks are collected in
        %   this internal buffer, after they are processed by the LDPC
        %   decoder core. The internal buffer may be reset using the reset
        %   method, when the retransmission of a particular information
        %   block is completed, allowing the transmission of another
        %   information block.
        b_hat_buffer
        
        %CODE_BLOCK_CRC_PASSED Vector of binary flags which indicates whether the
        %CRC check has been passed for each code block
        %   The vector of binary flags may be reset using the reset
        %   method, when the retransmission of a particular information
        %   block is completed, allowing the transmission of another
        %   information block.
        code_block_CRC_passed
        
    end
    
    methods
        % Constructor allowing properties to be set according to e.g.
        % a = NRLDPCDecoder('BG',1,'A',20,'G',132);
        function obj = NRLDPCDecoder(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Access = protected)
        
        function setupImpl(obj)
            B_ = obj.B;
            C_ = obj.C;
            N_cb_ = obj.N_cb;
            code_block_L_ = obj.code_block_L;
            
            obj.hTBCRCDetector = comm.CRCDetector('Polynomial',obj.transport_block_CRC_polynomial);
            if code_block_L_ > 0
                obj.hCBCRCDetector = comm.CRCDetector('Polynomial',obj.code_block_CRC_polynomial);
            end
%             try
%                 obj.hLDPCDecoder = comm.gpu.LDPCDecoder('ParityCheckMatrix',obj.H,'MaximumIterationCount',obj.iterations,'IterationTerminationCondition','Parity check satisfied');
%             catch
                obj.hLDPCDecoder = comm.LDPCDecoder('ParityCheckMatrix',obj.H,'MaximumIterationCount',obj.iterations,'IterationTerminationCondition','Parity check satisfied');
%             end
            d_tilde_buffer_ = cell(C_,1);
            for r=0:C_-1
                d_tilde_buffer_{r+1} = zeros(N_cb_,1);
            end
            obj.d_tilde_buffer = d_tilde_buffer_;
            
            obj.b_hat_buffer = zeros(B_,1);
            obj.code_block_CRC_passed = zeros(C_,1);
        end
        
        
        function a_hat = stepImpl(obj, g_tilde)
            f_tilde = code_block_concatenation(obj, g_tilde);
            e_tilde = bit_interleaving(obj, f_tilde);
            d_tilde = bit_selection(obj, e_tilde);
            c_hat = LDPC_coding(obj, d_tilde);
            b_hat = code_block_segmentation(obj, c_hat);
            a_hat = crc_calculation(obj, b_hat);
        end
        
        % Implements Section 5.5 of TS38.212
        function f_tilde = code_block_concatenation(obj, g_tilde)
            G_ = obj.G;
            C_ = obj.C;
            E_r_ = obj.E_r;
            
            if size(g_tilde,1) ~= G_ || size(g_tilde,2) ~= 1
                error('ldpc_3gpp_matlab:Error','g_tilde should be a column vector of length G.');
            end
            
            f_tilde = cell(C_,1);
            for r = 0:C_-1
                f_tilde{r+1} = zeros(E_r_(r+1), 1);
            end
            
            k = 0;
            r = 0;
            
            while r < C_
                j = 0;
                while j < E_r_(r+1)
                    f_tilde{r+1}(j+1) = g_tilde(k+1);
                    k = k + 1;
                    j = j + 1;
                end
                r = r + 1;
            end
        end
        
        % Implements Section 5.4.2.2 of TS38.212
        function e_tilde = bit_interleaving(obj, f_tilde)
            C_ = obj.C;
            E_r_ = obj.E_r;
            Q_m_ = obj.Q_m;
            
            if size(f_tilde,1) ~= C_ || size(f_tilde,2) ~= 1
                error('ldpc_3gpp_matlab:Error','f_tilde should be a column cell array of length C.');
            end
            
            e_tilde = cell(C_,1);
            
            for r = 0:C_-1
                
                if size(f_tilde{r+1},1) ~= E_r_(r+1) || size(f_tilde{r+1},2) ~= 1
                    error('ldpc_3gpp_matlab:Error','f_tilde{r+1} should be a column vector of length E_r(r+1).');
                end
                
                e_tilde{r+1} = zeros(E_r_(r+1),1);
                
                for j = 0:E_r_(r+1)/Q_m_-1
                    for i = 0:Q_m_-1
                        e_tilde{r+1}(i*E_r_(r+1)/Q_m_+j+1) = f_tilde{r+1}(i+j*Q_m_+1);
                    end
                end
            end
        end
        
        % Implements Section 5.4.2.1 of TS38.212
        function d_tilde = bit_selection(obj, e_tilde)
            C_ = obj.C;
            E_r_ = obj.E_r;
            N_ = obj.N;
            K_prime_ = obj.K_prime;
            Z_c_ = obj.Z_c;
            K_ = obj.K;
            k_0_ = obj.k_0;
            N_cb_ = obj.N_cb;
            I_HARQ_ = obj.I_HARQ;
            d_tilde_buffer_ = obj.d_tilde_buffer;
            
            if size(e_tilde,1) ~= C_ || size(e_tilde,2) ~= 1
                error('ldpc_3gpp_matlab:Error','e_tilde should be a column cell array of length C.');
            end
            
            d_tilde = cell(C_,1);
            
            for r=0:C_-1
                if size(e_tilde{r+1},1) ~= E_r_(r+1) || size(e_tilde{r+1},2) ~= 1
                    error('ldpc_3gpp_matlab:Error','e_tilde{r+1} should be a column vector of length E_r(r+1).');
                end
                
                d_tilde{r+1} = zeros(N_,1);
                d_tilde{r+1}(max(K_prime_-2*Z_c_+1,1):K_-2*Z_c_) = NaN;
                
                k = 0;
                j = 0;
                while k < E_r_(r+1)
                    if ~isnan(d_tilde{r+1}(mod(k_0_ + j, N_cb_)+1))
                        d_tilde{r+1}(mod(k_0_ + j, N_cb_)+1) = d_tilde{r+1}(mod(k_0_ + j, N_cb_)+1) + e_tilde{r+1}(k+1);
                        k = k+1;
                    end
                    j = j+1;
                end
                
                if I_HARQ_ ~= 0
                    d_tilde{r+1}(1:N_cb_) = d_tilde{r+1}(1:N_cb_) + d_tilde_buffer_{r+1};
                    d_tilde_buffer_{r+1} = d_tilde{r+1}(1:N_cb_);
                end
            end
            obj.d_tilde_buffer = d_tilde_buffer_;
        end
        
        % Implements Section 5.3.2 of TS38.212
        function c_hat = LDPC_coding(obj, d_tilde)
            C_ = obj.C;
            N_ = obj.N;
            Z_c_ = obj.Z_c;
            K_ = obj.K;
            
            if size(d_tilde,1) ~= C_ || size(d_tilde,2) ~= 1
                error('ldpc_3gpp_matlab:Error','d_tilde should be a column cell array of length C.');
            end
            
            c_hat = cell(C_,1);
            
            for r=0:C_-1
                if size(d_tilde{r+1},1) ~= N_ || size(d_tilde{r+1},2) ~= 1
                    error('ldpc_3gpp_matlab:Error','d_tilde{r+1} should be a column vector of length N.');
                end
                
                cw_tilde = [zeros(2*Z_c_,1); d_tilde{r+1}];
                c_tilde = cw_tilde(1:K_);
                cw_tilde(isnan(cw_tilde)) = inf;
                c_hat{r+1} = double(step(obj.hLDPCDecoder, cw_tilde));
                c_hat{r+1}(isnan(c_tilde)) = NaN;
            end
        end
        
        % Implements Section 5.2.2 of TS38.212
        function b_hat = code_block_segmentation(obj, c_hat)
            C_ = obj.C;
            I_HARQ_ = obj.I_HARQ;
            b_hat_buffer_ = obj.b_hat_buffer;
            B_ = obj.B;
            K_ = obj.K;
            K_prime_ = obj.K_prime;
            code_block_L_ = obj.code_block_L;
            CBGTI_flags_ = obj.CBGTI_flags;
            code_block_CRC_passed_ = obj.code_block_CRC_passed;
            
            if size(c_hat,1) ~= C_ || size(c_hat,2) ~= 1
                error('ldpc_3gpp_matlab:Error','c_hat should be a column cell array of length C.');
            end
            
            if I_HARQ_ ~= 0
                b_hat = b_hat_buffer_;
            else
                b_hat = zeros(B_,1);
            end
            
            s = 0;
            for r = 0:C_-1
                if size(c_hat{r+1},1) ~= K_ || size(c_hat,2) ~= 1
                    error('ldpc_3gpp_matlab:Error','c_hat should be a column vector of length K.');
                end
                
                code_block_CRC_failed = false;
                if C_ > 1
                    [~,code_block_CRC_failed] = step(obj.hCBCRCDetector,c_hat{r+1}(1:K_prime_));
                end
                
                for k = 0:K_prime_-code_block_L_-1
                    if ~code_block_CRC_failed && CBGTI_flags_(r+1) == 1
                        b_hat(s+1) = c_hat{r+1}(k+1);
                        code_block_CRC_passed_(r+1) = 1;
                    end
                    s = s + 1;
                end                
            end
            
            if I_HARQ_ ~= 0
                b_hat_buffer_ = b_hat;
            end
            
            obj.code_block_CRC_passed = code_block_CRC_passed_;
            obj.b_hat_buffer = b_hat_buffer_;
        end
        
        % Implements Section 5.1 of TS38.212
        function a_hat = crc_calculation(obj, b_hat)
            B_ = obj.B;
            A_ = obj.A;
            code_block_CRC_passed_ = obj.code_block_CRC_passed;
            
            if size(b_hat,1) ~= B_ || size(b_hat,2) ~= 1
                error('ldpc_3gpp_matlab:Error','b_hat should be a column vector of length B.');
            end
            
            a_hat = zeros(A_,1);
            
            for k = 0:A_-1
                a_hat(k+1) = b_hat(k+1);
            end
            
            [~,transport_block_CRC_failed] = step(obj.hTBCRCDetector,b_hat);
            if transport_block_CRC_failed || any(~code_block_CRC_passed_)
                a_hat = [];
            end
        end
        
        
        function resetImpl(obj)
            C_ = obj.C;
            B_ = obj.B;
            N_cb_ = obj.N_cb;
            
            d_tilde_buffer_ = cell(C_,1);
            for r=0:C_-1
                d_tilde_buffer_{r+1} = zeros(N_cb_,1);
            end
            obj.d_tilde_buffer = d_tilde_buffer_;
            
            obj.b_hat_buffer = zeros(B_,1);
            obj.code_block_CRC_passed = zeros(C_,1);            
        end
        
    end
end
