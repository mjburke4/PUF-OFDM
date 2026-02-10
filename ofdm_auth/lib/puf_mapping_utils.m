function out = puf_mapping_utils()
% Returns handles to small deterministic mapping utilities.
% Keep this file tiny and boring. Boring is reliable.

out.hash_bits_u32 = @hash_bits_u32;
out.xorshift32    = @xorshift32;
out.prg_bits      = @prg_bits;
out.mix_nonce     = @mix_challenge_nonce;
out.perm_from_seed= @perm_from_seed;
out.mask_from_seed= @mask_from_seed;

end

function s = hash_bits_u32(bits)
% bits: 0/1 row vector
% Simple FNV-like hash into uint32 (not crypto; just deterministic)
s = uint32(2166136261);
for i = 1:numel(bits)
    s = bitxor(s, uint32(bits(i)));
    s = uint32(mod(uint64(s) * 16777619, 2^32));
end
end

function s = xorshift32(s)
s = uint32(s);
s = bitxor(s, bitshift(s, 13));
s = bitxor(s, bitshift(s, -17));
s = bitxor(s, bitshift(s, 5));
end

function bits = prg_bits(seed_u32, nbits)
% Deterministic PRG -> 0/1 bits
s = uint32(seed_u32);
bits = false(1, nbits);
for i = 1:nbits
    s = xorshift32(s);
    bits(i) = bitand(s, uint32(1)) ~= 0;
end
bits = double(bits);
end

function Ceff = mix_challenge_nonce(Cseed, nonce_bits, target_len)
% Cseed, nonce_bits: 0/1 row vectors. Returns 0/1 length target_len.
% Deterministic mixing while preserving length.
utils = puf_mapping_utils();
seed = utils.hash_bits_u32([Cseed(:).' nonce_bits(:).']);
mix  = utils.prg_bits(seed, target_len);
Ceff = mod(Cseed(1:target_len) + mix, 2);
end

function p = perm_from_seed(seed_u32, N)
% Returns a permutation of 1..N using Fisher-Yates with xorshift PRNG
p = 1:N;
s = uint32(seed_u32);
for i = N:-1:2
    s = xorshift32(s);
    j = double(mod(s, uint32(i))) + 1;
    tmp = p(i); p(i) = p(j); p(j) = tmp;
end
end

function a = mask_from_seed(seed_u32, N, K)
% Binary mask length N with exactly K ones, deterministic
assert(K <= N);
p = perm_from_seed(seed_u32, N);
a = zeros(1, N);
a(p(1:K)) = 1;
end
