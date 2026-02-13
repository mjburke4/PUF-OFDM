function C_eff = mix_challenge_nonce(C_seed, nonce_bits, target_len)
% mix_challenge_nonce  Fixed-length challenge mixing with nonce freshness
% C_seed: 0/1 vector (any length)
% nonce_bits: 0/1 vector (any length)
% target_len: m

C_seed = logical(C_seed(:).');
nonce_bits = logical(nonce_bits(:).');

% Truncate/pad seed challenge
if numel(C_seed) < target_len
    C_seed = [C_seed false(1, target_len-numel(C_seed))];
else
    C_seed = C_seed(1:target_len);
end

if isempty(nonce_bits)
    C_eff = C_seed;
    return;
end

mask = prg_from_nonce(nonce_bits, target_len);
C_eff = xor(C_seed, mask);
end

function bits = prg_from_nonce(U, nbits)
% Simple deterministic PRG from nonce bits -> mask bits
% (Good enough for simulation; in real life you'd use a cryptographic KDF/PRF.)

s = uint32(2166136261);
U = double(U(:).');

for idx = 1:numel(U)
    s = bitxor(s, uint32(U(idx) + idx));
    s = uint32(mod(uint64(s) * uint64(16777619), 2^32));
end

bits = false(1, nbits);
x = s;
for j = 1:nbits
    x = bitxor(x, bitshift(x, 13, 'uint32'));
    x = bitxor(x, bitshift(x, -17, 'uint32'));
    x = bitxor(x, bitshift(x, 5, 'uint32'));
    bits(j) = bitand(x, uint32(1)) ~= 0;
    x = bitxor(x, uint32(2654435761));
end
end
