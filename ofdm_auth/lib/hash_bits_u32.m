function h = hash_bits_u32(bits)
% hash_bits_u32  Deterministic bits -> uint32 seed (FNV-ish style)
% bits: 0/1 vector

bits = double(bits(:).');
s = uint32(2166136261);

for k = 1:numel(bits)
    s = bitxor(s, uint32(bits(k) + k));
    s = uint32(mod(uint64(s) * uint64(16777619), 2^32));
end

h = s;
end
