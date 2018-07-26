classdef NRLDPC < matlab.System
    
    % Once the step function has been called, the values of nontunable
    % properties can only be changed if the release function is called
    % first.
    properties(Nontunable)
        
        % BG selects between LDPC base graph 1 of Table 5.3.2-2 in TS38.212 
        % and LDPC base graph 2 of Table 5.3.2-3. Base graph 1 is used for 
        % long block lengths and high coding rates, while base graph 2 is 
        % used for short block lengt hand low coding rates, as described in 
        % Sections 6.2.2 and 7.2.2 of TS38.212.
        BG = 1; % Default value
        
        % CRC selects between CRC16, CRC24A or CRC24B, as defined in 
        % Section 5.1 of TS38.212. Long transport blocks are appended with 
        % a CRC24A, while short transport blocks are appended with a CRC16, 
        % as described in Sections 6.2.1 and 7.2.1 of TS38.212. If the
        % transport block (with appended CRC) is sufficiently long, then it
        % is decomposed into two or more code blocks, each of which is 
        % appended with a CRC24B, as described in Section 5.2.2 of 
        % TS38.212. The use of CRC16 or CRC24 in this code implies that we 
        % are considering the LDPC coding of a transport block that is not 
        % long enough to be decomposed into two or more code blocks. The 
        % use of CRC24B in this code implies that we are considering one of 
        % the code blocks within a transport block that is long enough to 
        % be decomposed into two or more code blocks.
        CRC = 'CRC24B'; % Default value
                
        % K_prime_minus_L specifies the number of information bits. If we 
        % are considering the LDPC coding of a transport block that is not 
        % long enough to be decomposed into two or more code blocks, then 
        % the number of information bits in the transport block is given by 
        % A, as defined in Sections 6.2.1 and 6.3.1 of TS38.212. If we are 
        % considering one of the code blocks within a transport block that 
        % is long enough to be decomposed into two or more code blocks, 
        % then the number of information bits in the code block is given by 
        % K'-L, as defined in Section 5.2.2 of TS38.212.
        K_prime_minus_L = 20; % Default value
        
        % I_LBRM specifies whether a limit is imposed upon the lenghth of 
        % the circular buffer used for rate matching, as defined in Section 
        % 5.4.2.1 of TS38.212. A full buffer is used if I_LBRM = 0 and a 
        % limited buffer is used otherwise.
        I_LBRM = 0; % Default value
    end

    % Tunable properties can be changed anytime, even after the step 
    % function has been called.
    properties
        
        % N_ref specifies the limit imposed upon the lenghth of 
        % the circular buffer used for rate matching, when I_LBRM is
        % non-zero, as defined in Section 5.4.2.1 of TS38.212. N_ref is
        % ignored when I_LBRM is zero. N_ref is tunable so that it can be
        % changed for successive retransmissions during HARQ.
        N_ref = 132; % Default value

        % rv_id specifies the redundancy version number, as defined in
        % Section 5.4.2.1 of TS38.212. rv_id is tunable so that it can be
        % changed for successive retransmissions during HARQ.
        rv_id = 0; % Default value
        
        % E specifies the number of encoded bits in the output bit sequence
        % after rate matching, as defined in Section 5.4.2.1 of TS38.212. E
        % is tunable so that it can be changed for successive
        % retransmissions during HARQ.
        E = 132; % Default value
        
        % Q_m specifies the number of bits per PSK or QAM symbol, which is 
        % (against convention) referred to as 'modulation order' in Section
        % 5.4.2.2 of TS38.212. Q_m is tunable so that it can be changed for 
        % successive retransmissions during HARQ.
        Q_m = 1; % Default value
    end
        
    % Protected dependent properties cannot be set manually. Instead, they
    % are functions of the non-dependent properties.
    properties(Dependent, SetAccess = protected)
        
        % CRCPolynomial specifies the polynomial used when appending a CRC
        % to the information bits, as defined in Section 5.1 of TS38.212.
        CRCPolynomial
        
        % L specifies the length of the CRC appended to the information
        % bits, as defined in Section 5.1 of TS38.212. If we are 
        % considering the LDPC coding of a transport block that is not long 
        % enough to be decomposed into two or more code blocks, then the 
        % length of the CRC is given by L, as defined in Sections 6.2.1 and 
        % 6.3.1 of TS38.212. If we are considering one of the code blocks 
        % within a transport block that is long enough to be decomposed 
        % into two or more code blocks, then length of the CRC is given by 
        % L, as defined in Section 5.2.2 of TS38.212.
        L
        
        % K_prime specifies the total number of information and CRC bits,
        % as defined in Section 5.2.2 of TS38.212.
        K_prime
        
        % K_b specifies the value of the parameter K_b, as defined in
        % Section 5.2.2 of TS38.212.
        K_b
        
        % Z_c specifies the lifting size, as defined in Section 5.2.2 of
        % TS38.212.
        Z_c

        % K is the number of information, CRC and padding bits that are
        % LDPC coded, as defined in Section 5.2.2 of TS38.212.
        K

        % i_LS specifies the set index which contains Z_c, as defined in
        % Section 5.3.2 of TS38.212.
        i_LS
        
        % V is a matrix that specifies the rotations used by each non-zero 
        % element of the base graph specified by BG, for the case of the 
        % set index i_LS, as defined in Tables 5.3.2-2 and 5.3.2-3 of 
        % TS38.212. Note that the values in V are incremented by 1 relative 
        % to those in Tables 5.3.2-2 and  5.3.2-3 of TS38.212, since a 
        % sparse matrix cannot store the value 0 in Matlab.
        V
        
        % H is the parity check matrix that is obtained by lifting the base 
        % graph specified by BG using the lifting size Z_c and the 
        % rotations of V, as defined in Section 5.3.2 of TS38.212.
        H
        
        % N specifies the number of encoded bits, as defined in Section
        % 5.3.2 of TS38.212.
        N
        
        % N_cb specifies the length of the rate matching buffer, as
        % defined in Section 5.4.2.1 of TS38.212.
        N_cb
        
        % k_0 specifies the starting position of the redundancy version
        % rv_id, as defined in Table 5.4.2.1-2 of TS38.212.
        k_0
    end
    
    % A StringSet specifies the valid values of string parameters.
    properties (Hidden,Constant)
        CRCSet = matlab.system.StringSet({'CRC24A','CRC24B','CRC16'});
    end    
    
    % Methods used to set and get the values of properties. 
    methods
        
        % Constructor allowing properties to be set according to e.g.
        % a = NRLDPC('BG',1,'K_prime_minus_L',20,'E',132);
        function obj = NRLDPC(varargin)
            setProperties(obj,nargin,varargin{:});
        end
        
        
        function set.BG(obj, BG)
            if BG < 1 || BG > 2
                error('ldpc_3gpp_matlab:UnsupportedBaseGraph','Valid values of BG are 1 and 2.');
            end
            obj.BG = BG;
        end
        
        function set.K_prime_minus_L(obj, K_prime_minus_L)
            if K_prime_minus_L < 0
                error('ldpc_3gpp_matlab:UnsupportedBlockLength','K_prime_minus_L should not be negative.');
            end
            obj.K_prime_minus_L = K_prime_minus_L;
        end
        
        function set.N_ref(obj, N_ref)
            if N_ref < 0
                error('ldpc_3gpp_matlab:UnsupportedBlockLength','N_ref should not be negative.');
            end
            obj.N_ref = N_ref;
        end
        
        function set.rv_id(obj, rv_id)
            if rv_id < 0 || rv_id > 3
                error('ldpc_3gpp_matlab:UnsupportedParameter','Valid values of rv_id are 0, 1, 2 and 3.');
            end
            obj.rv_id = rv_id;
        end
        
        function set.E(obj, E)
            if E < 0
                error('ldpc_3gpp_matlab:UnsupportedBlockLength','E should not be negative.');
            end
            obj.E = E;
        end
        
        function set.Q_m(obj, Q_m)
            if Q_m <= 0
                error('ldpc_3gpp_matlab:UnsupportedBlockLength','Q_m should be positive.');
            end
            obj.Q_m = Q_m;
        end
        
        function CRCPolynomial = get.CRCPolynomial(obj)
            CRCPolynomial = get_3gpp_crc_polynomial(obj.CRC);
        end

        function L = get.L(obj)
            [~,L] = get_3gpp_crc_polynomial(obj.CRC);
        end

        function K_prime = get.K_prime(obj)
            K_prime = obj.K_prime_minus_L + obj.L;
        end
        
        function K_b = get.K_b(obj)
            if obj.BG == 1
                K_b = 22;
            elseif obj.BG == 2
                % TS 38.212 uses B rather than K_prime for the comparisons
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
                error('ldpc_3gpp_matlab:UnsupportedBaseGraph','BG must be 1 or 2');
            end
        end
        
        function Z_c = get.Z_c(obj)
            Z_c = get_3gpp_lifting_size(obj.K_b, obj.K_prime);
        end
        
        function i_LS = get.i_LS(obj)
            i_LS = get_3gpp_set_index(obj.Z_c);
        end
        
        function V = get.V(obj)
            V = get_3gpp_base_graph(obj.BG,obj.i_LS);
        end
        
        function H = get.H(obj)
            H = get_pcm(obj.V,obj.Z_c);
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
        
        function N_cb = get.N_cb(obj)
            if obj.I_LBRM == 0
                N_cb = obj.N;
            else
                N_cb = min(obj.N, obj.N_ref);
            end
        end
        
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
                    error('ldpc_3gpp_matlab:UnsupportedRedundancyVersion','Valid values of rv_id are 0, 1, 2 and 3.')
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
                    error('ldpc_3gpp_matlab:UnsupportedRedundancyVersion','Valid values of rv_id are 0, 1, 2 and 3.')
                end
            else
                error('ldpc_3gpp_matlab:UnsupportedBaseGraph','Valid values of BG are 1 and 2.');
            end
        end
        
        
    end
    
    
    methods(Access = protected)
        
        function setupImpl(obj)

        end
        
        function stepImpl(obj)

        end
       
        function resetImpl(obj)

        end
        
    end
end
