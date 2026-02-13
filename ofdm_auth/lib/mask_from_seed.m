function a = mask_from_seed(seed_u32, N, K)
% mask_from_seed  Deterministic K-of-N binary mask from uint32 seed

assert(K <= N, 'K must be <= N');

p = perm_from_seed(seed_u32, N);
a = zeros(1, N);
a(p(1:K)) = 1;
end
