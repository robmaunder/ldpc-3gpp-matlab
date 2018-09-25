%NRLDPC Base class for 3GPP New Radio LDPC objects
%   This base class cannot be used to perform the processing of any LDPC 
%   coding, but it can be used to setup all LDPC coding parameters. This 
%   may be inherited by derived classes that can perform the processing of 
%   LDPC coding.
%
%   Copyright © 2018 Robert G. Maunder. This program is free software: you 
%   can redistribute it and/or modify it under the terms of the GNU General 
%   Public License as published by the Free Software Foundation, either 
%   version 3 of the License, or (at your option) any later version. This 
%   program is distributed in the hope that it will be useful, but WITHOUT 
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
%   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
%   for more details.
classdef NRLDPC < matlab.System
    
    % Once the step function has been called, the values of nontunable
    % properties can only be changed if the release function is called
    % first.
    properties(Nontunable)
        
        %BG Basegraph selection
        %   Selects between LDPC base graph 1 of Table 5.3.2-2 in TS38.212 and LDPC 
        %   base graph 2 of Table 5.3.2-3. Base graph 1 is used for long block 
        %   lengths and high coding rates, while base graph 2 is used for short 
        %   block lengths and low coding rates, as described in Sections 6.2.2 and 
        %   7.2.2 of TS38.212.
        BG = 1; % Default value
        
        %CRC Cyclic Redundancy Check (CRC) selection
        %   Selects between 'CRC16', 'CRC24A' or 'CRC24B', as defined in Section
        %   5.1 of TS38.212. Alternatively, 'None' can be used to remove the CRC.
        %   Long transport blocks are appended with a CRC24A, while short transport
        %   blocks are appended with a CRC16, as described in Sections 6.2.1 and
        %   7.2.1 of TS38.212. If the transport block (with appended CRC) is
        %   sufficiently long, then it is decomposed into two or more code blocks,
        %   each of which is appended with a CRC24B, as described in Section 5.2.2
        %   of TS38.212. The use of CRC16 or CRC24 in this code implies that we are
        %   considering the LDPC coding of a transport block that is not long
        %   enough to be decomposed into two or more code blocks. The use of CRC24B
        %   in this code implies that we are considering one of the code blocks
        %   within a transport block that is long enough to be decomposed into two
        %   or more code blocks.
        CRC = 'CRC24B'; % Default value
                
        %K_PRIME_MINUS_L Number of information bits
        %   If we are considering the LDPC coding of a transport block that is not
        %   long enough to be decomposed into two or more code blocks, then the
        %   number of information bits in the transport block is given by A, as
        %   defined in Sections 6.2.1 and 6.3.1 of TS38.212. If we are considering
        %   one of the code blocks within a transport block that is long enough to 
        %   be decomposed into two or more code blocks, then the number of 
        %   information bits in the code block is given by K'-L, as defined in 
        %   Section 5.2.2 of TS38.212.
        K_prime_minus_L = 20; % Default value
        
        %I_LBRM Enable limited buffer rate matching
        %   Specifies whether or not a limit is imposed upon the lenghth of the
        %   circular buffer used for rate matching, as defined in Section 5.4.2.1
        %   of TS38.212. A full buffer is used if I_LBRM = 0 and a limited buffer
        %   is used otherwise.
        I_LBRM = 0; % Default value
    
        %N_REF Circular buffer limit
        %   Specifies limit imposed upon the lenghth of the circular buffer used
        %   for rate matching, when I_LBRM is non-zero, as defined in Section
        %   5.4.2.1 of TS38.212. N_ref is ignored when I_LBRM is zero.
        N_ref = 132; % Default value
    end

    % Tunable properties can be changed anytime, even after the step 
    % function has been called.
    properties
        
        %RV_ID Redundancy version number 
        %   Specifies the redundancy version number, as defined in Section 5.4.2.1
        %   of TS38.212. rv_id is tunable so that it can be changed for successive
        %   retransmissions during HARQ.
        rv_id = 0; % Default value
        
        %E Number of encoded bits        
        %   Specifies the number of encoded bits in the output bit sequence after
        %   rate matching, as defined in Section 5.4.2.1 of TS38.212. E is tunable
        %   so that it can be changed for successive retransmissions during HARQ.
        E = 132; % Default value
        
        %Q_M Number of bits per modulation symbol 
        %   Specifies the number of bits per PSK or QAM symbol, which is (against
        %   convention) referred to as 'modulation order' in Section 5.4.2.2 of
        %   TS38.212. Q_m is tunable so that it can be changed for successive
        %   retransmissions during HARQ.
        Q_m = 1; % Default value
    end
        
    % Protected dependent properties cannot be set manually. Instead, they
    % are calculated automatically as functions of the non-dependent 
    % properties.
    properties(Dependent, SetAccess = protected)
        
        %CRCPOLYNOMIAL Cyclic Redundancy Check (CRC) polynomial 
        %   Specifies the polynomial used when appending a CRC to the information
        %   bits, as defined in Section 5.1 of TS38.212.
        CRCPolynomial
        
        %L Cyclic Redundancy Check (CRC) length
        %   Specifies the length of the CRC appended to the information bits, as
        %   defined in Section 5.1 of TS38.212. If we are considering the LDPC
        %   coding of a transport block that is not long enough to be decomposed
        %   into two or more code blocks, then the length of the CRC is given by L,
        %   as defined in Sections 6.2.1 and 6.3.1 of TS38.212. If we are
        %   considering one of the code blocks within a transport block that is
        %   long enough to be decomposed into two or more code blocks, then length
        %   of the CRC is given by L, as defined in Section 5.2.2 of TS38.212.
        L
        
        %K_PRIME Number of information and CRC bits
        %   Specifies the total number of information and CRC bits, as defined in
        %   Section 5.2.2 of TS38.212.
        K_prime
        
        %K_B Parameter K_b
        %   Specifies the value of the parameter K_b, as defined in Section 5.2.2
        %   of TS38.212.
        K_b
        
        %Z_C Lifting size
        %   Specifies the lifting size, as defined in Section 5.2.2 of TS38.212.
        Z_c

        %K Number of information, CRC and padding bits
        %   Specifies the number of information, CRC and padding bits that are LDPC
        %   coded, as defined in Section 5.2.2 of TS38.212.
        K

        %I_LS Set index
        %   Specifies the set index which contains Z_c, as defined in Section 5.3.2
        %   of TS38.212.
        i_LS
        
        %V Rotations matrix        
        %   A matrix that specifies the rotations used by each non-zero element of
        %   the base graph specified by BG, for the case of the set index i_LS, as
        %   defined in Tables 5.3.2-2 and 5.3.2-3 of TS38.212. Note that the values
        %   in V are incremented by 1 relative to those in Tables 5.3.2-2 and
        %   5.3.2-3 of TS38.212, since a sparse matrix cannot store the value 0 in
        %   Matlab.
        V
        
        %H Parity Check Matrix (PCM)
        %   Specifies the PCM that is obtained by lifting the base graph specified
        %   by BG using the lifting size Z_c and the rotations of V, as defined in
        %   Section 5.3.2 of TS38.212.
        H
        
        %N Number of encoded bits        
        %   Specifies the number of encoded bits, as defined in Section 5.3.2 of
        %   TS38.212.
        N
        
        %N_CB Rate matching buffer length         
        %   Specifies the length of the rate matching buffer, as defined in Section
        %   5.4.2.1 of TS38.212.
        N_cb
        
        %K_0 Redundancy version starting position        
        %   Specifies the starting position of the redundancy version rv_id, as
        %   defined in Table 5.4.2.1-2 of TS38.212.
        k_0
    end
    
    % A StringSet specifies the valid values of string parameters.
    properties (Hidden,Constant)
        CRCSet = matlab.system.StringSet({'CRC24A','CRC24B','CRC16','None'});
    end    
    
    % Methods used to set and get the values of properties. 
    methods
        
        % Constructor allowing properties to be set according to e.g.
        % a = NRLDPC('BG',1,'K_prime_minus_L',20,'E',132);
        function obj = NRLDPC(varargin)
            setProperties(obj,nargin,varargin{:});
        end
        
        % Valid values of BG are described in Section 5.3.2 of TS38.212.
        function set.BG(obj, BG)
            if BG < 1 || BG > 2
                error('ldpc_3gpp_matlab:UnsupportedParameters','Valid values of BG are 1 and 2.');
            end
            obj.BG = BG;
        end
        
        function set.K_prime_minus_L(obj, K_prime_minus_L)
            if K_prime_minus_L < 0
                error('ldpc_3gpp_matlab:UnsupportedParameters','K_prime_minus_L should not be negative.');
            end
            obj.K_prime_minus_L = K_prime_minus_L;
        end
        
        function set.N_ref(obj, N_ref)
            if N_ref < 0
                error('ldpc_3gpp_matlab:UnsupportedParameters','N_ref should not be negative.');
            end
            obj.N_ref = N_ref;
        end
        
        % Valid values of rv_id are described in Section 5.4.2.1 of 
        % TS38.212.
        function set.rv_id(obj, rv_id)
            if rv_id < 0 || rv_id > 3
                error('ldpc_3gpp_matlab:UnsupportedParameters','Valid values of rv_id are 0, 1, 2 and 3.');
            end
            obj.rv_id = rv_id;
        end
        
        function set.E(obj, E)
            if E < 0
                error('ldpc_3gpp_matlab:UnsupportedParameters','E should not be negative.');
            end
            obj.E = E;
        end
        
        % Valid values of Q_m are derived from TS38.211.
        function set.Q_m(obj, Q_m)
            if Q_m ~= 1 && Q_m ~= 2 && Q_m ~= 4 && Q_m ~= 6 && Q_m ~= 8
                error('ldpc_3gpp_matlab:UnsupportedParameters','Valid vales of Q_m are 1, 2, 4, 6 and 8.');
            end
            obj.Q_m = Q_m;
        end
        
        % CRC polynomials are given in Section 5.1 of TS38.212.
        function CRCPolynomial = get.CRCPolynomial(obj)
            [CRCPolynomial,~] = get_3gpp_crc_polynomial(obj.CRC);
        end

        % CRC lengths are given in Section 5.1 of TS38.212.
        function L = get.L(obj)
            [~,L] = get_3gpp_crc_polynomial(obj.CRC);
        end

        function K_prime = get.K_prime(obj)
            K_prime = obj.K_prime_minus_L + obj.L;
        end
        
        % The calculation of K_b is given in Section 5.2.2 of TS38.212.
        function K_b = get.K_b(obj)
            if obj.BG == 1
                K_b = 22;
            elseif obj.BG == 2
                % TS38.212 uses B rather than K_prime for the comparisons
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
                error('ldpc_3gpp_matlab:UnsupportedParameters','BG must be 1 or 2');
            end
        end
        
        % The calculation of Z_c is given in Section 5.2.2 of TS38.212.
        function Z_c = get.Z_c(obj)
            Z_c = get_3gpp_lifting_size(obj.K_b, obj.K_prime);
        end
        
        % The calculation of K is given in Section 5.2.2 of TS38.212.
        function K = get.K(obj)
            if obj.BG == 1
                K = obj.Z_c*22;
            elseif obj.BG == 2
                K = obj.Z_c*10;
            else
                error('ldpc_3gpp_matlab:UnsupportedParameters','Valid values of BG are 1 and 2.');
            end
        end
        
        % The calculation of i_LS is given in Section 5.3.2 of TS38.212.
        function i_LS = get.i_LS(obj)
            i_LS = get_3gpp_set_index(obj.Z_c);
        end
        
        % The calculation of V is given in Section 5.3.2 of TS38.212.
        function V = get.V(obj)
            V = get_3gpp_base_graph(obj.BG,obj.i_LS);
        end
        
        % The calculation of H is given in Section 5.3.2 of TS38.212.
        function H = get.H(obj)
            H = get_pcm(obj.V,obj.Z_c);
        end
        
        % The calculation of N is given in Section 5.3.2 of TS38.212.
        function N = get.N(obj)
            if obj.BG == 1
                N = obj.Z_c*66;
            elseif obj.BG == 2
                N = obj.Z_c*50;
            else
                error('ldpc_3gpp_matlab:UnsupportedParameters','Valid values of BG are 1 and 2.');
            end
        end
        
        % The calculation of N_cb is given in Section 5.4.2.1 of TS38.212.
        function N_cb = get.N_cb(obj)
            if obj.I_LBRM == 0
                N_cb = obj.N;
            else
                N_cb = min(obj.N, obj.N_ref);
            end
        end
        
        % The calculation of k_0 is given in Table 5.4.2.1-2 of TS38.212.
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
                    error('ldpc_3gpp_matlab:UnsupportedParameters','Valid values of rv_id are 0, 1, 2 and 3.')
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
                    error('ldpc_3gpp_matlab:UnsupportedParameters','Valid values of rv_id are 0, 1, 2 and 3.')
                end
            else
                error('ldpc_3gpp_matlab:UnsupportedParameters','Valid values of BG are 1 and 2.');
            end
        end
    end
end
