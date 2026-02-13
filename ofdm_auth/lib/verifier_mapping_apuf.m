function map = verifier_mapping_apuf(node, C_seed, nonce_bits, B, Npilots, K_active)
% verifier_mapping_apuf  Derive per-symbol structure from APUF + nonce
%
% Inputs:
%   node      : node struct (distinct APUF instance)
%   C_seed    : base challenge bits (0/1)
%   nonce_bits: nonce bits (0/1)
%   B         : number of PUF response bits
%   Npilots   : number of pilot tones
%   K_active  : number of active pilot tones for IM
%
% Outputs:
%   map.C_eff
%   map.R_bits
%   map.seed_u32
%   map.pilot_perm_local   (1..Npilots)
%   map.pilot_mask_local   (0/1 length Npilots)

m = node.m;

% Effective challenge
C_eff = mix_challenge_nonce(C_seed, nonce_bits, m);

% PUF response bits (uses your uploaded get_puf_bits + arbiter_puf_sim_node)
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
