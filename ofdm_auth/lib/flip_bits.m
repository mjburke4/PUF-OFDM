function R2 = flip_bits(R, p)
% flip_bits  Flip each bit of R with probability p
% R: 0/1 vector
R = double(R(:).');
flips = rand(size(R)) < p;
R2 = xor(R, flips);
R2 = double(R2);
end