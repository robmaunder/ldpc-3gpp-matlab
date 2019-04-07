clear all

test_count = 0;
while true
    
    R = rand;
    I_LBRM = round(rand);
    TBS_LBRM = randi(100000);
    A = randi(100000);
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
    G = Q_m*N_L*round(randi(100000)/Q_m/N_L);
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
        encUL.LimitedBufferSize = TBS_LBRM;
    end
    setTransportBlock(encUL, a);
    Q_m_strings = {'pi/2-BPSK','QPSK','16QAM','64QAM','256QAM'};
    Q_m_string = Q_m_strings{Q_m_index};
    g_matlab = encUL(Q_m_string,N_L,G,rv_id);
    
    fprintf("test=%d\tBG=%d\tA=%d\tI_LBRM=%d\tTBS_LBRM=%d\trv_id=%d\tG=%d\tQ_m=%d\tN_L=%d\t", test_count,BG,A,I_LBRM,TBS_LBRM,rv_id,G,Q_m,N_L);
    if ~isequal(g,g_matlab)
        fprintf("ERROR\n");
    else
        fprintf("OK\n");
    end
    
    
    
    test_count = test_count+1;
end