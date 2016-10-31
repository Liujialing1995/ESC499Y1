function gene = encode_gene(phase)
% Maps a phase quantity to a 16-bit gray code 
    binary = uint16(mod(phase,2*pi)/(pi)*2^15);
    gene = binary;
    %gene = bitxor(binary, bitor(bitshift(binary,-1),bitshift(binary,31)));
end

