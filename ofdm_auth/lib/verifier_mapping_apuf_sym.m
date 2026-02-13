function map = verifier_mapping_apuf_sym(node, C_seed, nonce_frame, ell, B, Npilots, K_active)
% verifier_mapping_apuf_sym  Per-symbol mapping derived from APUF + (nonce_frame, ell)

m = node.m;

% Per-symbol nonce (frame nonce + ell bits)
nonce_bits = nonce_with_symbol(nonce_frame, ell, 16);

% Effective challenge for this symbol
C_eff = mix_challenge_nonce(C_seed, nonce_bits, m);

% PUF response bits (your CDMA-consistent mechanism)
R = get_puf_bits(node, C_eff, B);

% Seed
seed_u32 = hash_bits_u32(R);

% Mapping derived from seed
pilot_perm_local = perm_from_seed(seed_u32, Npilots);
pilot_mask_local = mask_from_seed(seed_u32, Npilots, K_active);

map.C_eff = C_eff;
map.R_bits = R;
map.seed_u32 = seed_u32;
map.pilot_perm_local = pilot_perm_local;
map.pilot_mask_local = pilot_mask_local;
end
