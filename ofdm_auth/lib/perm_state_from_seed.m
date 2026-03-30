function perm_auth = perm_state_from_seed(seed_u32)
% perm_state_from_seed
% Deterministic 2-state permutation on a 2-pilot auth subset.
%
% Output:
%   [1 2] -> identity
%   [2 1] -> swapped

    if mod(double(seed_u32), 2) == 0
        perm_auth = [1 2];
    else
        perm_auth = [2 1];
    end
end