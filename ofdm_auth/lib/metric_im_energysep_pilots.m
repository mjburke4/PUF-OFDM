function T = metric_im_energysep_pilots(Z, pilot_idx, Atilde)
% Energy separation only within pilot bins (pilot-domain IM)
% pilot_idx: vector of pilot bin indices
% Atilde: vector of predicted active pilot bin indices (subset of pilot_idx)

% ---- Defensive checks (catch drift / bugs early)
pilot_idx = pilot_idx(:);
Atilde    = Atilde(:);

assert(~isempty(pilot_idx), 'pilot_idx is empty.');
assert(~isempty(Atilde),    'Atilde is empty.');
assert(all(ismember(Atilde, pilot_idx)), ...
    'Atilde contains indices outside pilot_idx. Pilot-domain IM violated.');

assert(numel(pilot_idx) > numel(Atilde), ...
    'No inactive pilots remain: numel(pilot_idx) must be > numel(Atilde).');

e = abs(Z).^2;

pilot_mask  = false(size(Z));
pilot_mask(pilot_idx) = true;

active_mask = false(size(Z));
active_mask(Atilde) = true;

inactive_mask = pilot_mask & ~active_mask;

Ea = mean(e(active_mask));
Ei = mean(e(inactive_mask));

% ---- Raw statistic (your original)
T_raw = Ea - Ei;

% ---- Normalized statistic (recommended)
eps0 = 1e-12;
%T = (Ea - Ei) / (Ei + eps0);   % dimensionless, more SNR-robust
T = struct('raw', T_raw, 'norm', (Ea - Ei)/(Ei + eps0));

end
