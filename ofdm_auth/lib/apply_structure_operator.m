function [Xtilde, templ, st] = apply_structure_operator(cfg, X, Xp, seed_u32, meta)
% Applies permutation and/or pilot-domain IM mask in frequency domain.
% Returns:
%   Xtilde: physical symbol (data+pilots after structure)
%   templ : physical pilot template for Metric A
%   st    : struct describing active set, permutation, etc.

utils = puf_mapping_utils();

N = cfg.Nfft;
pilot_idx = meta.pilot_idx(:).';

% --- Pilot-domain IM mask (only on pilot tones)
aP = ones(1, numel(pilot_idx));
if isfield(cfg,'K_active_pilots') && cfg.K_active_pilots < numel(pilot_idx)
    aP = utils.mask_from_seed(seed_u32, numel(pilot_idx), cfg.K_active_pilots);
end

% --- Pilot permutation (within pilot set only)
perm_depth = min(cfg.perm_depth, numel(pilot_idx));
p_local = 1:numel(pilot_idx);
p_draw  = utils.perm_from_seed(seed_u32, numel(pilot_idx));
p_local(1:perm_depth) = p_draw(1:perm_depth);

% Apply: permute pilot positions among pilot bins, then mask
Xtilde = X;

% Extract logical pilots as a vector
pilot_vals = Xp(pilot_idx);

% Permute pilots within pilot bins
pilot_vals_perm = pilot_vals(p_local);

% Mask (inactive pilots -> 0)
pilot_vals_perm = pilot_vals_perm .* aP(:);

% Write back into pilot bins
Xtilde(pilot_idx) = pilot_vals_perm;

% Template for correlator is pilot-only physical
templ = zeros(N,1);
templ(pilot_idx) = pilot_vals_perm;

% Active physical set for IM metric (within pilots)
active_local = find(aP > 0);
Atilde = pilot_idx(active_local);

st.aP = aP;
st.p_local = p_local;
st.Atilde = Atilde(:);
st.pilot_idx = pilot_idx(:);
end
