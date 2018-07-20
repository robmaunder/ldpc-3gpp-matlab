clear all

valid_lifting_sizes = get_3gpp_valid_lifting_sizes;

modulation_order = 4;
hMod = comm.PSKModulator(modulation_order, 'BitInput',true);


while true
    BG = randi([1 2]);
    if BG == 1
        K_prime = randi([4,8448]);
    else
        K_prime = randi([4,3840]);
    end


    E = floor(K_prime/rand());

    
    rv_id = 0; %randi([0,3]);

    I_LBRM = randi([0,1]);

    
    hEnc = NRLDPCEncoder('BG',BG,'K_prime',K_prime,'I_LBRM',I_LBRM,'E',E,'rv_id',rv_id);
    hDec = NRLDPCDecoder('BG',BG,'K_prime',K_prime,'I_LBRM',I_LBRM,'E',E,'rv_id',rv_id,'iterations',10);
    hDec2 = NRLDPCDecoder2('BG',BG,'K_prime',K_prime,'I_LBRM',I_LBRM,'E',E,'rv_id',rv_id,'iterations',10);
    
    N = hEnc.N;

    N_ref = randi([K_prime,N]);
    hEnc.N_ref = N_ref;
    hDec.N_ref = N_ref;
    hDec2.N_ref = N_ref;

    


    EsN0 = 20;
    
    fprintf("%d\t%d\t%d\t%d\t%d\t%d\t%d\n", BG, K_prime, I_LBRM, hEnc.N_cb, E, rv_id, EsN0)
    
    
    
    hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',EsN0);
    hDemod = comm.PSKDemodulator(modulation_order, 'BitOutput',true, 'DecisionMethod','Log-likelihood ratio', 'Variance', 1/10^(hChan.SNR/10));
    
    
    b = round(rand(K_prime,1));
    c = [b;NaN(hEnc.K-K_prime,1)];
    d = step(hEnc,c);
    e = d(~isnan(d));
    f = [e;zeros(mod(-length(e),log2(modulation_order)),1)];
    tx = step(hMod, f);
    rx = step(hChan, tx);
    f_tilde = step(hDemod, rx);
    e_tilde = f_tilde(1:length(e));
    d_tilde = d;
    d_tilde(~isnan(d)) = e_tilde;
    c_hat = step(hDec, d_tilde);
    b_hat = c_hat(~isnan(c_hat));
    c_hat2 = step(hDec2, d_tilde);
    b_hat2 = c_hat2(~isnan(c_hat2));
    

    if ~isequal(b_hat, b_hat2)
        fprintf("Decoder 1 mismatched with Decoder 2\n");
    end
    
    if ~isequal(b_hat, b)
        fprintf("Decoder 1 mismatched with Encoder\n");
    end
    
    if ~isequal(b_hat2, b)
        fprintf("Decoder 2 mismatched with Encoder\n");
    end
    
    
end





