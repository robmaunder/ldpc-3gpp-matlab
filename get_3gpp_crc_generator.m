function [crc_generator, L] = get_3gpp_crc_generator(crc_generator_id)

    if strcmp(crc_generator_id, 'CRC24A')
        crc_generator = 'z^24 + z^23 + z^18 + z^17 + z^14 + z^11 + z^10 + z^7 + z^6 + z^5 + z^4 + z^3 + z + 1';
        L = 24;
    elseif strcmp(crc_generator_id, 'CRC24B')
        crc_generator = 'z^24 + z^23 + z^6 + z^5 + z + 1';
        L = 24;
    elseif strcmp(crc_generator_id, 'CRC16')
        crc_generator = 'z^16 + z^12 + z^5 + 1';
        L = 16;
    else
        error('ldpc_3gpp_matlab:UnsupportedCRCGenerator','Invalid CRC generator identifier.');
    end


end
