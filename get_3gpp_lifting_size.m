function Z_c = get_3gpp_lifting_size(K_b, K_prime)

    valid_lifting_sizes = get_3gpp_valid_lifting_sizes();
    
    Z_c = inf;

    for set_index = 0:length(valid_lifting_sizes)-1
        Z_c_set = min(valid_lifting_sizes{set_index+1}(K_b*valid_lifting_sizes{set_index+1} >= K_prime));
        if Z_c_set < Z_c
            Z_c = Z_c_set;
        end
    end
    
    if Z_c == inf
        error('ldpc_3gpp_matlab:UnsupportedParameters','Invalid block length.');
    end
end
        