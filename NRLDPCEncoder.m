%NRLDPCENCODER 3GPP New Radio LDPC encoder
classdef NRLDPCEncoder < NRLDPC
    
    properties(Access = private, Hidden)
        %HTBCRCGENERATOR Cyclic Redundancy Check (CRC) generator
        %   A COMM.CRCGENERATOR used to generate the transport block CRC.  
        %
        %   See also COMM.CRCGENERATOR
        hTBCRCGenerator

        %HCBCRCGENERATOR Cyclic Redundancy Check (CRC) generator
        %   A COMM.CRCGENERATOR used to generate the code block CRC.  
        %
        %   See also COMM.CRCGENERATOR
        hCBCRCGenerator

        %HLDPCENCODER Low Density Parity Check (LDPC) encoder
        %   A COMM.LDPCENCODER used to perform the LDPC encoding.  
        %
        %   See also COMM.LDPCENCODER        
        hLDPCEncoder
    end
    
    methods
        % Constructor allowing properties to be set according to e.g.
        % a = NRLDPCEncoder('BG',1,'A',20,'G',132);
        function obj = NRLDPCEncoder(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end    
    
    % Methods used to execute processing.
    methods(Access = protected)
        
        % Code executed on the first time that the step function is called,
        % or the first time after the release function is called. e.g.
        % a = NRLDPCEncoder;
        % step(a); % <- setupImpl executed here
        % step(a); % <- setupImpl not executed here
        % reset(a);
        % step(a); % <- setupImpl not executed here
        % release(a);
        % step(a); % <- setupImpl executed here
        function setupImpl(obj)
            obj.hTBCRCGenerator = comm.CRCGenerator('Polynomial',obj.transport_block_CRC_polynomial);
            if obj.code_block_L > 0
                obj.hCBCRCGenerator = comm.CRCGenerator('Polynomial',obj.code_block_CRC_polynomial);
            end
            obj.hLDPCEncoder = comm.LDPCEncoder('ParityCheckMatrix',obj.H);
        end
        
        % Code executed by the step function. e.g.
        % a = NRLDPCEncoder;
        % step(a); % <- stepImpl executed here
        % step(a); % <- stepImpl executed here
        % reset(a);
        % step(a); % <- stepImpl executed here
        % release(a);
        % step(a); % <- stepImpl executed here
        function g = stepImpl(obj, a)
            b = crc_calculation(obj, a);
            c = code_block_segmentation(obj, b);
            d = LDPC_coding(obj, c);
            e = bit_selection(obj, d);
            f = bit_interleaving(obj, e);
            g = code_block_concatenation(obj, f);
        end
        
        % Implements Section 5.3.2 of TS38.212
        function b = crc_calculation(obj, a)
            if size(a,1) ~= obj.A || size(a,2) ~= 1
                error('ldpc_3gpp_matlab:Error','a should be a column vector of length A.');
            end
            
            b = zeros(obj.B,1);
            
            ap = step(obj.hTBCRCGenerator, a);
            p = ap(obj.A+1:obj.B);
                
            for k = 0:obj.A-1
                b(k+1) = a(k+1);
            end
            for k = obj.A:obj.B-1
                b(k+1) = p(k-obj.A+1);
            end
        end
                
        % Implements Section 5.2.2 of TS38.212
        function c = code_block_segmentation(obj, b)
            if size(b,1) ~= obj.B || size(b,2) ~= 1
                error('ldpc_3gpp_matlab:Error','b should be a column vector of length B.');
            end
            
            c = cell(obj.C,1);
            
            s = 0;           
            for r = 0:obj.C-1
                c{r+1} = zeros(obj.K,1);
                
                for k = 0:obj.K_prime-obj.code_block_L-1
                    c{r+1}(k+1) = b(s+1);
                    s = s + 1;
                end
                if obj.C > 1
                    cp = step(obj.hCBCRCGenerator, c{r+1}(1:obj.K_prime-obj.code_block_L));
                    p = cp(obj.K_prime-obj.code_block_L+1:obj.K_prime);
                    for k = obj.K_prime-obj.code_block_L:obj.K_prime-1
                        c{r+1}(k+1) = p(k+obj.code_block_L-obj.K_prime+1);
                    end
                end
                for k = obj.K_prime:obj.K-1
                    c{r+1}(k+1) = NaN;
                end
            end
        end
        
        % Implements Section 5.3.2 of TS38.212
        function d = LDPC_coding(obj, c)
            if size(c,1) ~= obj.C || size(c,2) ~= 1
                error('ldpc_3gpp_matlab:Error','c should be a column cell array of length C.');
            end
            
            d = cell(obj.C,1);
            
            for r = 0:obj.C-1
                if size(c{r+1},1) ~= obj.K || size(c{r+1},2) ~= 1
                    error('ldpc_3gpp_matlab:Error','c{r+1} should be a column vector of length K.');
                end
                
                d{r+1} = zeros(obj.N,1);
                
                % Not sure about what to do if there are NaNs within the first
                % 2*obj.Z_c elements of c. The following code (adapted from
                % TS38.212) does not set these to 0.
                for k = 2*obj.Z_c:obj.K-1
                    if ~isnan(c{r+1}(k+1))
                        d{r+1}(k-2*obj.Z_c+1) = c{r+1}(k+1);
                    else
                        c{r+1}(k+1) = 0;
                        d{r+1}(k-2*obj.Z_c+1) = NaN;
                    end
                end

                cw = step(obj.hLDPCEncoder, c{r+1});
                w = cw(obj.K+1:obj.N+2*obj.Z_c);

                for k = obj.K:obj.N+2*obj.Z_c-1
                    d{r+1}(k-2*obj.Z_c+1) = w(k-obj.K+1);
                end
            end
        end
        
        % Implements Section 5.4.2.1 of TS38.212
        function e = bit_selection(obj, d)
            if size(d,1) ~= obj.C || size(d,2) ~= 1
                error('ldpc_3gpp_matlab:Error','d should be a column cell array of length C.');
            end
            
            e=cell(obj.C,1);
            
            for r = 0:obj.C-1
                if size(d{r+1},1) ~= obj.N || size(d{r+1},2) ~= 1
                    error('ldpc_3gpp_matlab:Error','d{r+1} should be a column vector of length N.');
                end
                e{r+1} = zeros(obj.E_r(r+1),1);

                k = 0;
                j = 0;
                while k < obj.E_r(r+1)
                    if ~isnan(d{r+1}(mod(obj.k_0 + j, obj.N_cb)+1))
                        e{r+1}(k+1) = d{r+1}(mod(obj.k_0 + j, obj.N_cb)+1);
                        k = k+1;
                    end
                    j = j+1;
                end
            end
        end
        
        % Implements Section 5.4.2.2 of TS38.212
        function f = bit_interleaving(obj, e)
            if size(e,1) ~= obj.C || size(e,2) ~= 1
                error('ldpc_3gpp_matlab:Error','e should be a column cell array of length C.');
            end
            
            f = cell(obj.C,1);
            
            for r = 0:obj.C-1
                
                if size(e{r+1},1) ~= obj.E_r(r+1) || size(e{r+1},2) ~= 1
                    error('ldpc_3gpp_matlab:Error','e{r+1} should be a column vector of length E_r(r+1).');
                end
            
                f{r+1} = zeros(obj.E_r(r+1),1);

                for j = 0:obj.E_r(r+1)/obj.Q_m-1
                    for i = 0:obj.Q_m-1
                        f{r+1}(i+j*obj.Q_m+1) = e{r+1}(i*obj.E_r(r+1)/obj.Q_m+j+1);
                    end
                end
            end
        end

        % Implements Section 5.5 of TS38.212
        function g = code_block_concatenation(obj, f)
            if size(f,1) ~= obj.C || size(f,2) ~= 1
                error('ldpc_3gpp_matlab:Error','f should be a column cell array of length C.');
            end
            for r = 0:obj.C-1                
                if size(f{r+1},1) ~= obj.E_r(r+1) || size(f{r+1},2) ~= 1
                    error('ldpc_3gpp_matlab:Error','f{r+1} should be a column vector of length E_r(r+1).');
                end
            end
            
            g = zeros(obj.G,1);
            
            k = 0;
            r = 0;
            
            while r < obj.C
                j = 0;
                while j < obj.E_r(r+1)
                    g(k+1) = f{r+1}(j+1);
                    k = k + 1;
                    j = j + 1;
                end
                r = r + 1;
            end
        end
        
        % Code executed by the reset function. e.g.
        % a = NRLDPCEncoder;
        % step(a);
        % step(a);
        % reset(a); % <- resetImpl executed here
        % step(a);
        % release(a);
        % step(a);
        function resetImpl(obj)

        end
        
    end
end
