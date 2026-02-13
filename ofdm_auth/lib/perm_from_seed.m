function p = perm_from_seed(seed_u32, N)
% p = perm_from_seed(seed_u32, N)
%
% Deterministic permutation of 1..N derived from a uint32 seed.
% Uses a local xorshift32 PRNG (does NOT touch MATLAB global RNG).
%
% Inputs:
%   seed_u32 : uint32 seed (e.g., from hash_bits_u32(PUF_bits))
%   N        : permutation size
%
% Output:
%   p        : 1xN permutation of 1..N

if ~isa(seed_u32, 'uint32')
    seed_u32 = uint32(seed_u32);
end

% Initialize permutation vector
p = 1:N;

% Local PRNG state
s = seed_u32;
if s == 0
    s = uint32(1);  % avoid zero-lock state
end

% Fisher-Yates shuffle (descending)
for i = N:-1:2
    s = xorshift32_local(s);

    % Draw j in [1, i]
    j = double(mod(s, uint32(i))) + 1;

    % Swap p(i) <-> p(j)
    tmp = p(i);
    p(i) = p(j);
    p(j) = tmp;
end

end


% ---------------- Local PRNG ----------------
function s = xorshift32_local(s)
% One round of xorshift32
s = bitxor(s, bitshift(s, 13, 'uint32'));
s = bitxor(s, bitshift(s, -17, 'uint32'));
s = bitxor(s, bitshift(s, 5, 'uint32'));
end
