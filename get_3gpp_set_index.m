function i_LS = get_3gpp_set_index(Z)

% Table 5.3.2-1 of TS38.212
valid_lifting_sizes = {
    [2, 4, 8, 16, 32, 64, 128, 256]
    [3, 6, 12, 24, 48, 96, 192, 384]
    [5, 10, 20, 40, 80, 160, 320]
    [7, 14, 28, 56, 112, 224]
    [9, 18, 36, 72, 144, 288]
    [11, 22, 44, 88, 176, 352]
    [13, 26, 52, 104, 208]
    [15, 30, 60, 120, 240]
    };

i_LS = 0;
while i_LS < 8 && isempty(find(valid_lifting_sizes{i_LS+1}==Z, 1))
    i_LS = i_LS+1;
end
if i_LS == 8
    error('ldpc_3gpp_matlab:UnsupportedLifingSize','Invalid lifting size.');
end

end

