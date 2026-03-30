function a = mask_from_seed_subset(seed_u32, N, K)
% Deterministic K-of-N mask on a subset of pilots

    assert(K <= N, 'K must be <= N');

    p = perm_from_seed(seed_u32, N);
    a = zeros(1, N);
    a(p(1:K)) = 1;
end