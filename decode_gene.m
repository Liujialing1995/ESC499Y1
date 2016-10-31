function phase = decode_gene(gc)
% Maps a 16-bit gray code to a phase quantity in the range of -pi to pi.
%     gc = bitxor(gc, bitor(bitshift(gc,-16),bitshift(gc,16)));
%     gc = bitxor(gc, bitor(bitshift(gc,-8),bitshift(gc,24)));
%     gc = bitxor(gc, bitor(bitshift(gc,-4),bitshift(gc,28)));
%     gc = bitxor(gc, bitor(bitshift(gc,-2),bitshift(gc,30)));
%     gc = bitxor(gc, bitor(bitshift(gc,-1),bitshift(gc,31)));
    phase = double(gc)/double(2^15) * pi;
    if(phase > pi)
        phase = phase - 2*pi;
    end
end
