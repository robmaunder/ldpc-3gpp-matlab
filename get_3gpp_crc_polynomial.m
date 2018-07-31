function [crc_polynomial, L] = get_3gpp_crc_polynomial(crc)

    if strcmp(crc, 'CRC24A')
        crc_polynomial = 'z^24 + z^23 + z^18 + z^17 + z^14 + z^11 + z^10 + z^7 + z^6 + z^5 + z^4 + z^3 + z + 1';
        L = 24;
    elseif strcmp(crc, 'CRC24B')
        crc_polynomial = 'z^24 + z^23 + z^6 + z^5 + z + 1';
        L = 24;
    elseif strcmp(crc, 'CRC16')
        crc_polynomial = 'z^16 + z^12 + z^5 + 1';
        L = 16;
    elseif strcmp(crc, 'None')
        crc_polynomial = '';
        L=0;
    else
        error('ldpc_3gpp_matlab:UnsupportedParameters','Invalid CRC identifier.');
    end

end