classdef NRLDPCEncoder < NRLDPC
    
    properties(Access = private, Hidden)
        hCRCGenerator
        hLDPCEncoder
    end
    
    methods
        function obj = NRLDPCEncoder(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end    
    
    methods(Access = protected)
        
        function setupImpl(obj)
            obj.hCRCGenerator = comm.CRCGenerator('Polynomial',obj.CRCPolynomial);
            obj.hLDPCEncoder = comm.LDPCEncoder('ParityCheckMatrix',obj.ParityCheckMatrix);
        end
        
        function f = stepImpl(obj, b)
            c = append_CRC_and_padding(obj, b);
            d = LDPC_coding(obj, c);
            e = bit_selection(obj, d);
            f = bit_interleaving(obj, e);
        end
        
        function c = append_CRC_and_padding(obj, b)
            if length(b) ~= obj.K_prime_minus_L
                error('ldpc_3gpp_matlab:Error','Length of b should be K_prime_minus_L.');
            end
            
            c = zeros(obj.K, 1);
            
            s = 0;
            for k = 0:obj.K_prime_minus_L-1
                c(k+1) = b(s+1);
                s = s + 1;
            end
            bp = step(obj.hCRCGenerator, b);
            p = bp(obj.K_prime_minus_L+1:obj.K_prime);
            for k = obj.K_prime_minus_L:obj.K_prime-1
                c(k+1) = p(k+obj.L-obj.K_prime+1);
            end
            for k = obj.K_prime:obj.K-1
                c(k+1) = NaN;
            end
        end
        
        % Implements Section 5.3.2 of TS38.212
        function d = LDPC_coding(obj, c)
            if length(c) ~= obj.K
                error('ldpc_3gpp_matlab:Error','Length of c should be K.');
            end
            
            d = zeros(obj.N,1);
            
            % Not sure about what to do if there are NaNs within the first
            % 2*obj.Z_c elements of c. The following code (adapted from
            % TS38.212) does not set these to 0.
            for k = 2*obj.Z_c:obj.K-1
                if ~isnan(c(k+1))
                    d(k-2*obj.Z_c+1) = c(k+1);
                else
                    c(k+1) = 0;
                    d(k-2*obj.Z_c+1) = NaN;
                end
            end
            
            cw = step(obj.hLDPCEncoder, c);
            w = cw(obj.K+1:obj.N+2*obj.Z_c);
            
            for k = obj.K:obj.N+2*obj.Z_c-1
                d(k-2*obj.Z_c+1) = w(k-obj.K+1);
            end
        end
        
        % Implements Section 5.4.2.1 of TS38.212
        function e = bit_selection(obj, d)
            if length(d) ~= obj.N
                error('ldpc_3gpp_matlab:Error','Length of d should be N.');
            end
            
            e = zeros(obj.E,1);
            
            k = 0;
            j = 0;
            while k < obj.E
                if ~isnan(d(mod(obj.k_0 + j, obj.N_cb)+1))
                    e(k+1) = d(mod(obj.k_0 + j, obj.N_cb)+1);
                    k = k+1;
                end
                j = j+1;
            end
        end
        
        % Implements Section 5.4.2.2 of TS38.212
        function f = bit_interleaving(obj, e)
            if length(e) ~= obj.E
                error('ldpc_3gpp_matlab:Error','Length of e should be E.');
            end
            
            f = zeros(obj.E,1);
            
            for j = 0:obj.E/obj.Q_m-1
                for i = 0:obj.Q_m-1
                    f(i+j*obj.Q_m+1) = e(i*obj.E/obj.Q_m+j+1);
                end
            end
        end
        
        
        function resetImpl(obj)

        end
        
    end
end
