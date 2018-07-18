function H = get_pcm(H_BG, Z)

H = spalloc(size(H_BG,1)*Z, size(H_BG,2)*Z, nnz(H_BG)*Z);

[rows,cols,circ] = find(H_BG);

for edge_index = 1:length(circ)
    H(rows(edge_index)*Z+(1-Z:0), cols(edge_index)*Z+(1-Z:0)) = circshift(speye(Z),mod(circ(edge_index)-1,Z),2);
end

end

