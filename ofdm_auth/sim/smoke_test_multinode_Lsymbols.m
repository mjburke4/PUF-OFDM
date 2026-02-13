function smoke_test_multinode_Lsymbols()
cfg = cfg_default();

% Defaults requested
U = 10; m = 64; B = 64;

% Use L from cfg (or override here)
L = cfg.Lsym;

% Node population compatible with arbiter_puf_sim_node (needs field A)
nodes = init_node_population(U, m, 'seed', 7);   % your updated init_node_population creates nd.A

% Base challenge and per-frame nonce
C_seed = randi([0 1], 1, m);
nonce_frame = randi([0 1], 1, 32);

% Choose transmitting node
u0 = 3;

% Thresholds (pick something reasonable; tune later)
tauP  = 0.7;      % for |Tperm|
tauIM = 0.5;      % for normalized IM (depends on your definition)

% Accumulators per hypothesis user
Sp = zeros(1,U);   % soft sum of |T|
Sim   = zeros(1,U);   % soft sum of Tim_norm

Vperm = zeros(1,U);   % vote counts
Vim   = zeros(1,U);

% Channel config for this test
SNRdB = 10;
cfg.chan.type = 'awgn';

for ell = 1:L
    % Logical OFDM symbol
    [X, Xp, meta] = gen_ofdm_symbol(cfg);

    % --- TX mapping for this symbol (depends on ell)
    map_tx = verifier_mapping_apuf_sym(nodes(u0), C_seed, nonce_frame, ell, B, cfg.Npilots, cfg.K_active_pilots);
    [Xtilde, ~, ~] = apply_structure_operator_pilots(X, Xp, meta, map_tx.pilot_perm_local, map_tx.pilot_mask_local);

    % OFDM modulate
    x = ifft(Xtilde, cfg.Nfft);
    x_cp = [x(end-cfg.Ncp+1:end); x];

    % Channel + AWGN
    [y_cp, H] = channel_model(cfg, x_cp);
    sig_pow = mean(abs(y_cp).^2);
    noise_pow = sig_pow / (10^(SNRdB/10));
    w = sqrt(noise_pow/2) * (randn(size(y_cp)) + 1j*randn(size(y_cp)));
    r_cp = y_cp + w;

    % RX FFT + perfect equalization
    r = r_cp(cfg.Ncp+1:end);
    Y = fft(r, cfg.Nfft);
    Hhat = estimate_channel(cfg, H, [], [], meta.pilot_idx);
    Z = Y ./ (Hhat + 1e-12);

    % Test all user hypotheses
    for u = 1:U
        map_rx = verifier_mapping_apuf_sym(nodes(u), C_seed, nonce_frame, ell, B, cfg.Npilots, cfg.K_active_pilots);
        [~, templ_rx, st_rx] = apply_structure_operator_pilots(X, Xp, meta, map_rx.pilot_perm_local, map_rx.pilot_mask_local);

        % Metric A
        tP = abs(metric_perm_correlator(Z, templ_rx));
        Sp(u) = Sp(u) + tP;
        Vperm(u) = Vperm(u) + (tP >= tauP);

        % Metric B (pilot-only IM)
        tIM = metric_im_energysep_pilots(Z, meta.pilot_idx, st_rx.Atilde);
        % If your metric returns struct, use tIM = tIM.norm;
        if isstruct(tIM), tIM = tIM.norm; end

        Sim(u) = Sim(u) + tIM;
        Vim(u) = Vim(u) + (tIM >= tauIM);
    end
end

% Winners by soft combining
[~, best_perm_soft] = max(Sp);
[~, best_im_soft]   = max(Sim);

% Winners by hard voting
[~, best_perm_vote] = max(Vperm);
[~, best_im_vote]   = max(Vim);

fprintf('TX user u0=%d, SNR=%g dB, L=%d symbols\n', u0, SNRdB, L);

T = table((1:U).', Sp.', Vperm.', Sim.', Vim.', ...
    'VariableNames', {'User','Sp_sum','Vperm_votes','Sim_sum','Vim_votes'});
disp(T);

fprintf('Best (perm soft) user: %d\n', best_perm_soft);
fprintf('Best (perm vote) user: %d\n', best_perm_vote);
fprintf('Best (IM soft)   user: %d\n', best_im_soft);
fprintf('Best (IM vote)   user: %d\n', best_im_vote);
end
