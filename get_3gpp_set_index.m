function i_LS = get_3gpp_set_index(Z)

valid_lifting_sizes = get_3gpp_valid_lifting_sizes();

i_LS = 0;
while i_LS < 8 && isempty(find(valid_lifting_sizes{i_LS+1}==Z, 1))
    i_LS = i_LS+1;
end
if i_LS == 8
    error('ldpc_3gpp_matlab:UnsupportedParameters','Invalid lifting size.');
end

end

