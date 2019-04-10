%   This script compares NRLDPCEncoder with nrULSCH from the 5G Toolbox and
%   confirms that both are compliant with each other and with 3GPP TS
%   38.212.
%
%
%   Copyright © 2018 Robert G. Maunder. This program is free software: you
%   can redistribute it and/or modify it under the terms of the GNU General
%   Public License as published by the Free Software Foundation, either
%   version 3 of the License, or (at your option) any later version. This
%   program is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
%   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
%   for more details.


clear all

test_count = 0;
while true
    
    R = rand;
    I_LBRM = round(rand);
    A = ceil(100000^rand);
    TBS_LBRM = round(A/rand);
    a = round(rand(A,1));
    if A<=292 || (A <= 3824 && R <= 0.67) || R<=0.25
        BG = 2;
    else
        BG = 1;
    end
    Q_m_values = [1 2 4 6 8];
    Q_m_index = randi(5);
    Q_m = Q_m_values(Q_m_index);
    N_L = randi(4);
    G = Q_m*N_L*round(A/rand/Q_m/N_L);
    rv_id = randi([0 3]);
    
    
    enc = NRLDPCEncoder;
    enc.BG = BG;
    enc.A = A;
    enc.I_LBRM = I_LBRM;
    enc.TBS_LBRM = TBS_LBRM;
    enc.rv_id = rv_id;
    enc.G = G;
    enc.Q_m = Q_m;
    enc.N_L = N_L;
    try
        g = enc(a);
    catch ME
        if strcmp(ME.identifier, 'ldpc_3gpp_matlab:UnsupportedParameters')
            continue
        else
            rethrow(ME);
        end
    end
    
    
    encUL = nrULSCH;
    encUL.TargetCodeRate = R;
    encUL.LimitedBufferRateMatching = I_LBRM == 1;
    if I_LBRM == 1
        encUL.LimitedBufferSize = enc.N_ref;
    end
    setTransportBlock(encUL, a);
    Q_m_strings = {'pi/2-BPSK','QPSK','16QAM','64QAM','256QAM'};
    Q_m_string = Q_m_strings{Q_m_index};
    g_matlab = encUL(Q_m_string,N_L,G,rv_id);
    
    fprintf("test=%d\tBG=%d\tA=%d\tI_LBRM=%d\tTBS_LBRM=%d\trv_id=%d\tG=%d\tQ_m=%d\tN_L=%d\n", test_count,BG,A,I_LBRM,TBS_LBRM,rv_id,G,Q_m,N_L);
    if ~isequal(g,g_matlab)
        error('Mismatch!');
    end
    
    
    
    test_count = test_count+1;
end