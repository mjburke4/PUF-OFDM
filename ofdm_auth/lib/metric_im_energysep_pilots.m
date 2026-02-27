function T = metric_im_energysep_pilots(Z, pilot_idx, Atilde)
% Pilot-domain IM energy separation metric with robust normalization.
% Z: N x 1 equalized (or raw FFT) bins
% pilot_idx: pilot bin indices (vector)
% Atilde: predicted active pilot indices (subset of pilot_idx)

pilot_idx = pilot_idx(:);
Atilde = Atilde(:);

assert(all(ismember(Atilde, pilot_idx)), 'Atilde must be subset of pilot_idx.');
assert(numel(pilot_idx) > numel(Atilde), 'Need at least one inactive pilot.');

e = abs(Z).^2;

pilot_mask = false(size(Z));
pilot_mask(pilot_idx) = true;

active_mask = false(size(Z));
active_mask(Atilde) = true;

inactive_mask = pilot_mask & ~active_mask;

Ea = mean(e(active_mask));
Ei = mean(e(inactive_mask));

eps0 = 1e-12;
T = (Ea - Ei) / (Ea + Ei + eps0);   % in roughly [-1,1]
end
