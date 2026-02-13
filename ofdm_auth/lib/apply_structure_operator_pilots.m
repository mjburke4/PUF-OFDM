function [Xtilde, templ, st] = apply_structure_operator_pilots(X, Xp, meta, pilot_perm_local, pilot_mask_local)
% Apply pilot permutation + pilot-domain IM mask using provided perm/mask.
% X: full logical OFDM spectrum
% Xp: pilot-only logical vector
% meta.pilot_idx: pilot subcarrier indices (length Npilots)

pilot_idx = meta.pilot_idx(:).';
Npilots = numel(pilot_idx);

assert(numel(pilot_perm_local) == Npilots);
assert(numel(pilot_mask_local) == Npilots);

pilot_vals = Xp(pilot_idx);

% permute pilots within pilot bins
pilot_vals_perm = pilot_vals(pilot_perm_local);

% apply mask
pilot_vals_perm = pilot_vals_perm .* pilot_mask_local(:);

% write back
Xtilde = X;
Xtilde(pilot_idx) = pilot_vals_perm;

templ = zeros(size(X));
templ(pilot_idx) = pilot_vals_perm;

st.Atilde = pilot_idx(pilot_mask_local(:) > 0);
st.pilot_idx = pilot_idx(:);
st.pilot_perm_local = pilot_perm_local(:);
st.pilot_mask_local = pilot_mask_local(:);
end
