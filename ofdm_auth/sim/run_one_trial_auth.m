function out = run_one_trial_auth(cfg, nodes, u0, SNRref_dB, is_legit, p_flip)
if nargin < 6, p_flip = 0; end
% function out = run_one_trial_auth(cfg, nodes, u0, SNRref_dB, is_legit, p_flip)
% if nargin < 6, p_flip = 0; end
%function out = run_one_trial_auth(cfg, nodes, u0, SNRref_dB, is_legit)
% Returns accept decisions for both metrics for one frame
L = cfg.Lsym;
m = nodes(u0).m;
B = 64;
d0 = 10; gamma = 2.7;

C_seed = randi([0 1], 1, m);
nonce_frame = randi([0 1], 1, 32);

u_tx = u0;   % default legitimate TX identity % TX identity (used only for pathloss distance)

u_tx = u0;   % default legitimate transmitter

if ~is_legit
    switch cfg.adv_mode
        case 'impostor_othernode'
            cand = setdiff(1:numel(nodes), u0);
            assert(~isempty(cand), 'Need at least 2 nodes for impostor_othernode.');
            u_tx = cand(randi(numel(cand)));

        case 'impostor_random'
            u_tx = u0;   % claimed identity is still u0, but structure will be random

        case 'replay_stale_nonce'
            u_tx = u0;   % same node, stale nonce mismatch

        otherwise
            error('Unknown cfg.adv_mode: %s', cfg.adv_mode);
    end
end

% U = numel(nodes);
% L = cfg.Lsym;
% m = nodes(u0).m;
% B = 64;   % your default
% d0 = 10; gamma = 2.7;

% Nonce/challenge per frame
% C_seed = randi([0 1], 1, m);
% nonce_frame = randi([0 1], 1, 32);

%u_tx = u0;  % TX identity (used only for pathloss distance)

% Noise floor calibrated at d0 using actual signal power each symbol
Vperm = 0;
Vim = 0;

Tp = zeros(1,L);
Tim = zeros(1,L);

% fprintf('DEBUG: is_legit=%d, adv_mode=%s, u0=%d, u_tx=%d\n', ...
%     is_legit, cfg.adv_mode, u0, u_tx);

for ell = 1:L
    [X, Xp, meta] = gen_ofdm_symbol(cfg);

    % --- TX-side effective challenge / nonce ---
    nonce_bits_tx = nonce_with_symbol(nonce_frame, ell, 16);
    C_eff_tx      = mix_challenge_nonce(C_seed, nonce_bits_tx, nodes(u0).m);

    if is_legit
        % Legitimate transmitter uses enrolled node u0
        R_tx = get_puf_bits(nodes(u0), C_eff_tx, B);
        R_tx = flip_bits(R_tx, p_flip);

        seed_tx = hash_bits_u32(R_tx);
        pilot_perm_local = perm_from_seed(seed_tx, cfg.Npilots);
        pilot_mask_local = mask_from_seed(seed_tx, cfg.Npilots, cfg.K_active_pilots);

    else
        switch cfg.adv_mode
            case 'impostor_random'
                % Random structure impostor
                seed_tx = uint32(randi([0, 2^32-1]));
                pilot_perm_local = perm_from_seed(seed_tx, cfg.Npilots);
                pilot_mask_local = mask_from_seed(seed_tx, cfg.Npilots, cfg.K_active_pilots);

            case 'impostor_othernode'
                % Different enrolled node pretending to be u0
                R_tx = get_puf_bits(nodes(u_tx), C_eff_tx, B);
                R_tx = flip_bits(R_tx, p_flip);

                seed_tx = hash_bits_u32(R_tx);
                pilot_perm_local = perm_from_seed(seed_tx, cfg.Npilots);
                pilot_mask_local = mask_from_seed(seed_tx, cfg.Npilots, cfg.K_active_pilots);

            case 'replay_stale_nonce'
                % Same node, but stale/wrong nonce at TX
                nonce_frame_tx = [1 - nonce_frame(1), nonce_frame(2:end)];
                nonce_bits_tx  = nonce_with_symbol(nonce_frame_tx, ell, 16);
                C_eff_tx       = mix_challenge_nonce(C_seed, nonce_bits_tx, nodes(u0).m);

                R_tx = get_puf_bits(nodes(u0), C_eff_tx, B);
                R_tx = flip_bits(R_tx, p_flip);

                seed_tx = hash_bits_u32(R_tx);
                pilot_perm_local = perm_from_seed(seed_tx, cfg.Npilots);
                pilot_mask_local = mask_from_seed(seed_tx, cfg.Npilots, cfg.K_active_pilots);

            otherwise
                error('Unknown cfg.adv_mode: %s', cfg.adv_mode);
        end
    end

    [Xtilde, ~, ~] = apply_structure_operator_pilots( ...
        X, Xp, meta, pilot_perm_local, pilot_mask_local);
    
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
if isempty(cfg.tauP)
    error('cfg.tauP is empty — must be set before running trial');
end

if isempty(cfg.tauIM)
    error('cfg.tauIM is empty — must be set before running trial');
end

out.accept_perm = frame_decision_votes(Vperm, L, cfg.alpha);
out.accept_im   = frame_decision_votes(Vim,   L, cfg.alpha);
out.accept_hybrid = out.accept_perm && out.accept_im;
out.u_tx = u_tx;

out.Tp = Tp;
out.Tim = Tim;

end
