function nonce_sym = nonce_with_symbol(nonce_frame, ell, ell_bits)
% nonce_with_symbol  Append symbol index bits to a frame nonce.
% nonce_frame : 0/1 row/col vector
% ell         : symbol index (1..L)
% ell_bits    : number of bits for ell encoding (default 16)

if nargin < 3, ell_bits = 16; end

nf = logical(nonce_frame(:).');
ell_u = uint32(ell);
ell_vec = de2bi(ell_u, ell_bits, 'left-msb');   % 0/1 row
nonce_sym = [nf ell_vec];
end
