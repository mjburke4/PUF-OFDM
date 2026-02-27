function out = run_one_trial_auth(cfg, nodes, u0, SNRref_dB, is_legit, p_flip)
if nargin < 6, p_flip = 0; end
%function out = run_one_trial_auth(cfg, nodes, u0, SNRref_dB, is_legit)
% Returns accept decisions for both metrics for one frame

U = numel(nodes);
L = cfg.Lsym;
m = nodes(u0).m;
B = 64;   % your default
d0 = 10; gamma = 2.7;

% Nonce/challenge per frame
C_seed = randi([0 1], 1, m);
nonce_frame = randi([0 1], 1, 32);

u_tx = u0;  % TX identity (used only for pathloss distance)

% Noise floor calibrated at d0 using actual signal power each symbol
Vperm = 0;
Vim = 0;

Tp = zeros(1,L);
Tim = zeros(1,L);

for ell = 1:L
    [X, Xp, meta] = gen_ofdm_symbol(cfg);
        % --- TX mapping
    if is_legit
        % Legit TX uses the enrolled mapping for u0 (must match verifier hypothesis)
        %map_tx = verifier_mapping_apuf_sym(nodes(u0), C_seed, nonce_frame, ell, B, cfg.Npilots, cfg.K_active_pilots);
        %[Xtilde, ~, ~] = apply_structure_operator_pilots(X, Xp, meta, map_tx.pilot_perm_local, map_tx.pilot_mask_local);
        % --- TX computes effective challenge
        C_eff = mix_challenge_nonce(C_seed, nonce_with_symbol(nonce_frame, ell, 16), nodes(u0).m);
        
        % --- TX reads PUF bits, then suffers flips
        R = get_puf_bits(nodes(u0), C_eff, B);
        R_noisy = flip_bits(R, p_flip);
        
        % --- TX derives seed and mapping from noisy bits
        seed_u32 = hash_bits_u32(R_noisy);
        pilot_perm_local = perm_from_seed(seed_u32, cfg.Npilots);
        pilot_mask_local = mask_from_seed(seed_u32, cfg.Npilots, cfg.K_active_pilots);
        
        [Xtilde, ~, ~] = apply_structure_operator_pilots(X, Xp, meta, pilot_perm_local, pilot_mask_local);
    else
        % Impostor TX uses random mapping per symbol (no access to PUF mapping)
        seed_u32 = uint32(randi([0, 2^32-1]));
        pilot_perm_local = perm_from_seed(seed_u32, cfg.Npilots);
        pilot_mask_local = mask_from_seed(seed_u32, cfg.Npilots, cfg.K_active_pilots);
        [Xtilde, ~, ~] = apply_structure_operator_pilots(X, Xp, meta, pilot_perm_local, pilot_mask_local);
    end

    
    x = ifft(Xtilde, cfg.Nfft);
    x_cp = [x(end-cfg.Ncp+1:end); x];

    % Calibrate noise at reference distance d0
    Pref = mean(abs(x_cp).^2);
    noise_pow = Pref / (10^(SNRref_dB/10));

    % Pathloss scaling for TX user
    amp = sqrt((d0 / nodes(u_tx).distance_m)^gamma);
    xchan_in = amp * x_cp;

    [y_cp, ~] = channel_model(cfg, xchan_in);

    w = sqrt(noise_pow/2) * (randn(size(y_cp)) + 1j*randn(size(y_cp)));
    r_cp = y_cp + w;

    r = r_cp(cfg.Ncp+1:end);
    Y = fft(r, cfg.Nfft);

    % No equalization (AWGN only), keep amplitude effects in observation
    Z = Y;

    % Verifier hypothesis ALWAYS tests u0
    map_rx = verifier_mapping_apuf_sym(nodes(u0), C_seed, nonce_frame, ell, B, cfg.Npilots, cfg.K_active_pilots);
    [~, templ_rx, st_rx] = apply_structure_operator_pilots(X, Xp, meta, map_rx.pilot_perm_local, map_rx.pilot_mask_local);
    % if is_legit && ell==1
    % fprintf('LEGIT seed tx=%u rx=%u\n', map_tx.seed_u32, map_rx.seed_u32);
    % end

    % Metric A: active-pilot correlator (best version)
    %tP = abs(metric_perm_correlator_activepilots(Z, templ_rx, st_rx.Atilde));
    Tc = metric_perm_correlator_activepilots(Z, templ_rx, st_rx.Atilde); % complex
    tP = real(Tc);   % coherence-sensitive


    % Metric B: normalized pilot-domain IM (Ea-Ei)/(Ea+Ei)
    tIM = metric_im_energysep_pilots(Z, meta.pilot_idx, st_rx.Atilde);
    if isstruct(tIM), tIM = tIM.norm; end

    % Per-symbol vote
    Vperm = Vperm + (tP >= cfg.tauP);
    Vim   = Vim   + (tIM >= cfg.tauIM);

    Tp(ell) = tP;
    Tim(ell) = tIM;

end

out.Vperm = Vperm;
out.Vim = Vim;

% if ~is_legit
%     fprintf('Impostor votes this frame: %d\n', Vperm);
% end

out.accept_perm = frame_decision_votes(Vperm, L, cfg.alpha);
out.accept_im   = frame_decision_votes(Vim,   L, cfg.alpha);
out.accept_hybrid = out.accept_perm && out.accept_im;
out.u_tx = u_tx;

out.Tp = Tp;
out.Tim = Tim;

%fprintf('permutation mean: %d std: 5d\n', mean(out.Tp), std(out.Tp));
%fprintf('IM mean: %d std: 5d\n', mean(out.Tim), std(out.Tim));
end
