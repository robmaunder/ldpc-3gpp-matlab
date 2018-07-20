valid_lifting_sizes = get_3gpp_valid_lifting_sizes;

hMod = comm.PSKModulator(2, 'BitInput',true);


while true
    BG = randi([1 2]);
    set_index = randi([0,7]);
    Z_c = valid_lifting_sizes{set_index+1}(randi([1,length(valid_lifting_sizes{set_index+1})]));
    
    hEnc = NRLDPCEncoder('BG',BG,'Z_c',Z_c);
    hDec = NRLDPCDecoder('BG',BG,'Z_c',Z_c,'iterations',10);
    hDec2 = NRLDPCDecoder2('BG',BG,'Z_c',Z_c,'iterations',10);
    K = hEnc.K;
    K_prime = randi([2*Z_c,K]);
    hEnc.K_prime = K_prime;
    hDec.K_prime = K_prime;
    hDec2.K_prime = K_prime;
    
    EsN0 = randi([-10, 10]);
    
    fprintf("%d\t%d\t%d\t%d\t\n", BG, Z_c, K_prime, EsN0)
    
    
    
    hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',EsN0);
    hDemod = comm.PSKDemodulator(2, 'BitOutput',true, 'DecisionMethod','Log-likelihood ratio', 'Variance', 1/10^(hChan.SNR/10));
    
    
    b = round(rand(K_prime,1));
    c = [b;NaN(K-K_prime,1)];
    d = step(hEnc,c);
    e = d(~isnan(d));
    
    tx = step(hMod, e);
    rx = step(hChan, tx);
    e_tilde = step(hDemod, rx);
    d_tilde = d;
    d_tilde(~isnan(d)) = e_tilde;
    c_hat = step(hDec, d_tilde);
    b_hat = c_hat(~isnan(c_hat));
    c_hat2 = step(hDec2, d_tilde);
    b_hat2 = c_hat2(~isnan(c_hat2));
    
    if ~isequal(b_hat, b_hat2)
        fprintf("Error\n");
        
        [b_hat, b_hat2]
    end
    
    
    
    
end





