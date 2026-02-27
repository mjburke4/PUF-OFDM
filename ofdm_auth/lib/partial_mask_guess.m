function a_guess = partial_mask_guess(a_true, rho)
% partial_mask_guess  Match a fraction rho of 1-entries in a_true, randomize remaining
% a_true: 0/1 length N
% rho in [0,1]

a_true = a_true(:).';
N = numel(a_true);
idx1 = find(a_true==1);
K = numel(idx1);

Kmatch = round(rho*K);
match_idx = idx1(randperm(K, Kmatch));

a_guess = zeros(1,N);
a_guess(match_idx) = 1;

% Fill remaining ones randomly among the zeros not already chosen
remain = K - Kmatch;
avail = setdiff(find(a_true==0), match_idx);
if remain > 0
    pick = avail(randperm(numel(avail), remain));
    a_guess(pick) = 1;
end
end
