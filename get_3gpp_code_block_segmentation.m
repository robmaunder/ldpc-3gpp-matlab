% Implements the first half of Section 5.2.2 of TS38.212

function [C,L,K_prime,K,Z_c] = get_3gpp_code_block_segmentation(B, BG)
if B <= 0
    error('ldpc_3gpp_matlab:UnsupportedBlockLength','B should be greater than 0.');
end

if BG == 1
    K_cb = 8448;
    K_b = 22;
elseif BG == 2
    K_cb = 3840;
    if B > 640
        K_b = 10;
    elseif B > 560
        K_b = 9;
    elseif B > 192
        K_b = 8;
    else
        K_b = 6;
    end
else
    error('ldpc_3gpp_matlab:UnsupportedBaseGraph','BG must be 1 or 2');
end

if B <= K_cb
    L = 0;
    C = 1;
    B_prime = B;
else
    L = 24;
    C = ceil(B/(K_cb-L));
    B_prime = B + C*L;
end

K_prime = B_prime/C;

% Not sure what to do if K_prime is not an integer

Z_c = get_3gpp_lifting_size(K_b, K_prime);

if BG == 1
    K = 22*Z_c;
elseif BG == 2
    K = 10*Z_c;
else
    error('ldpc_3gpp_matlab:UnsupportedBaseGraph','BG must be 1 or 2');
end
