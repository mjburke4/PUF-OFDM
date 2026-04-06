function fig_replay_far_vs_snr_all_clean_pretty_v2()

%% =========================
% Common settings
% ==========================
alpha = 0.75;
SNR_list = [-5 5 15 25];
MC = 6000; %2000;

% Calibration settings
SNR_cal   = 15;
L_cal     = 16; %8;
MC_cal    = 3000; %3000;
targetFAR = 0.01;
p_cal     = 0;

U = 10;
m = 64;
u0 = 3;

nodes = init_node_population(U, m, 'seed', 7);

% Fix distance so SNRref is meaningful for perm/IM branch
d0 = 10;
for u = 1:U
    nodes(u).distance_m = d0;
end

%% =========================
% Perm / IM config
% ==========================
cfg = cfg_default();
cfg.chan.type = 'awgn';
cfg.alpha = alpha;
cfg.Lsym  = L_cal;

tauP_grid  = linspace(0.00, 0.50, 101);
tauIM_grid = linspace(0.00, 1.00, 121);

% --- Calibrate tauP
best_tauP = tauP_grid(end);
best_FRR_P = inf;
best_FAR_P = inf;

for tau = tauP_grid
    cfgTmp = cfg;
    cfgTmp.tauP = tau;

    acc = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, false, 0);
        acc = acc + o.accept_perm;
    end
    FAR_here = acc / MC_cal;

    rej = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, true, p_cal);
        rej = rej + (~o.accept_perm);
    end
    FRR_here = rej / MC_cal;

    if (FAR_here <= targetFAR) && (FRR_here < best_FRR_P)
        best_tauP = tau;
        best_FRR_P = FRR_here;
        best_FAR_P = FAR_here;
    end
end

% --- Calibrate tauIM
best_tauIM = tauIM_grid(end);
best_FRR_IM = inf;
best_FAR_IM = inf;

for tau = tauIM_grid
    cfgTmp = cfg;
    cfgTmp.tauIM = tau;

    acc = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, false, 0);
        acc = acc + o.accept_im;
    end
    FAR_here = acc / MC_cal;

    rej = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, true, p_cal);
        rej = rej + (~o.accept_im);
    end
    FRR_here = rej / MC_cal;

    if (FAR_here <= targetFAR) && (FRR_here < best_FRR_IM)
        best_tauIM = tau;
        best_FRR_IM = FRR_here;
        best_FAR_IM = FAR_here;
    end
end

cfg.tauP  = best_tauP;
cfg.tauIM = best_tauIM;

%% =========================
% DPPA config
% ==========================
cfgD = cfg_default_dppa();
cfgD.chan.type = 'awgn';
cfgD.alpha = alpha;
cfgD.Lsym  = L_cal;

% Optional WLAN consistency:
% cfgD.Npilots = 4;

% Strongly consider using richer phase alphabet for replay robustness:
cfgD.phi_choices_deg = [-28 -10 10 28];

tauD_grid = linspace(0.60, 1.00, 81);

best_tauD = tauD_grid(end);
best_FRR_D = inf;
best_FAR_D = inf;

for tau = tauD_grid
    cfgTmp = cfgD;
    cfgTmp.tauD = tau;

    acc = 0;
    for t = 1:MC_cal
        o = run_one_trial_dppa(cfgTmp, nodes, u0, SNR_cal, false, 0);
        acc = acc + o.accept_dppa;
    end
    FAR_here = acc / MC_cal;

    rej = 0;
    for t = 1:MC_cal
        o = run_one_trial_dppa(cfgTmp, nodes, u0, SNR_cal, true, p_cal);
        rej = rej + (~o.accept_dppa);
    end
    FRR_here = rej / MC_cal;

    if (FAR_here <= targetFAR) && (FRR_here < best_FRR_D)
        best_tauD = tau;
        best_FRR_D = FRR_here;
        best_FAR_D = FAR_here;
    end
end

cfgD.tauD = best_tauD;

fprintf('\nReplay FAR figure calibration @ SNR=%g dB, L=%d:\n', SNR_cal, L_cal);
fprintf('  tauP  = %.4f   (FAR=%.4f, FRR=%.4f)\n', cfg.tauP,  best_FAR_P,  best_FRR_P);
fprintf('  tauIM = %.4f   (FAR=%.4f, FRR=%.4f)\n', cfg.tauIM, best_FAR_IM, best_FRR_IM);
fprintf('  tauD  = %.4f   (FAR=%.4f, FRR=%.4f)\n\n', cfgD.tauD, best_FAR_D, best_FRR_D);

%% =========================
% Replay evaluation
% ==========================
replay_FAR_P   = zeros(size(SNR_list));
replay_FAR_IM  = zeros(size(SNR_list));
replay_FAR_D   = zeros(size(SNR_list));
replay_FAR_HYB = zeros(size(SNR_list));

for si = 1:numel(SNR_list)
    snr = SNR_list(si);

    accP = 0; accIM = 0; accD = 0; accHYB = 0;

    for t = 1:MC
        nonce_old = randi([0 1], 1, 32);
        nonce_new = randi([0 1], 1, 32);

        outPI = run_one_trial_replay_pi(cfg, nodes, u0, snr, nonce_old, nonce_new);
        outD  = run_one_trial_replay_dppa(cfgD, nodes, u0, snr, nonce_old, nonce_new);

        accP   = accP   + outPI.accept_perm;
        accIM  = accIM  + outPI.accept_im;
        accD   = accD   + outD.accept_dppa;
        accHYB = accHYB + (outPI.accept_perm && outPI.accept_im);
    end

    replay_FAR_P(si)   = accP   / MC;
    replay_FAR_IM(si)  = accIM  / MC;
    replay_FAR_D(si)   = accD   / MC;
    replay_FAR_HYB(si) = accHYB / MC;

    fprintf('SNR=%g dB: Replay FAR(P)=%.4f FAR(IM)=%.4f FAR(D)=%.4f FAR(HYB)=%.4f\n', ...
        snr, replay_FAR_P(si), replay_FAR_IM(si), replay_FAR_D(si), replay_FAR_HYB(si));
end

%% =========================
% Pretty plot
% ==========================
figure('Color','w', 'Position',[120 120 900 580]);
hold on; grid on; box on;

set(gca, 'YScale','log', ...
         'FontSize',13, ...
         'LineWidth',1.1, ...
         'XMinorGrid','off', ...
         'YMinorGrid','on');

plot(SNR_list, max(replay_FAR_P,   1/MC), '-o', 'LineWidth',2.2, 'MarkerSize',8);
plot(SNR_list, max(replay_FAR_IM,  1/MC), '-s', 'LineWidth',2.2, 'MarkerSize',8);
plot(SNR_list, max(replay_FAR_D,   1/MC), '-d', 'LineWidth',2.2, 'MarkerSize',8);
plot(SNR_list, max(replay_FAR_HYB, 1/MC), '-^', 'LineWidth',2.2, 'MarkerSize',8);

xlabel('SNR_{ref} (dB)', 'FontSize', 14);
ylabel('Probability (log scale)', 'FontSize', 14);
title(sprintf('Replay Robustness: Perm vs IM vs DPPA vs Hybrid (L = %d, \\alpha = %.2f)', ...
    cfg.Lsym, cfg.alpha), ...
    'FontSize', 15, 'FontWeight', 'bold');

legend('Replay FAR (Perm)', 'Replay FAR (IM)', 'Replay FAR (DPPA)', 'Replay FAR (Hybrid)', ...
    'Location','northeast', 'FontSize', 12);

xticks(SNR_list);
xlim([min(SNR_list) max(SNR_list)]);
ylim([1/MC 1]);

end

% ==========================================================
function out = run_one_trial_replay_pi(cfg, nodes, u0, SNRref_dB, nonce_old, nonce_new)

L = cfg.Lsym;
m = nodes(u0).m;
B = 64;

C_seed = randi([0 1], 1, m);

Vperm = 0;
Vim   = 0;

for ell = 1:L
    [X, Xp, meta] = gen_ofdm_symbol(cfg);

    map_tx = verifier_mapping_apuf_sym(nodes(u0), C_seed, nonce_old, ell, B, ...
        cfg.Npilots, cfg.K_active_pilots);
    [Xtilde, ~, ~] = apply_structure_operator_pilots( ...
        X, Xp, meta, map_tx.pilot_perm_local, map_tx.pilot_mask_local);

    x = ifft(Xtilde, cfg.Nfft);
    x_cp = [x(end-cfg.Ncp+1:end); x];

    Pref = mean(abs(x_cp).^2);
    noise_pow = Pref / (10^(SNRref_dB/10));

    [y_cp, ~] = channel_model(cfg, x_cp);
    w = sqrt(noise_pow/2) * (randn(size(y_cp)) + 1j*randn(size(y_cp)));
    r_cp = y_cp + w;

    r = r_cp(cfg.Ncp+1:end);
    Y = fft(r, cfg.Nfft);
    Z = Y;

    map_rx = verifier_mapping_apuf_sym(nodes(u0), C_seed, nonce_new, ell, B, ...
        cfg.Npilots, cfg.K_active_pilots);
    [~, templ_rx, st_rx] = apply_structure_operator_pilots( ...
        X, Xp, meta, map_rx.pilot_perm_local, map_rx.pilot_mask_local);

    Tc = metric_perm_correlator_activepilots(Z, templ_rx, st_rx.Atilde);
    tP = real(Tc);

    tIM = metric_im_energysep_pilots(Z, meta.pilot_idx, st_rx.Atilde);

    Vperm = Vperm + (tP  >= cfg.tauP);
    Vim   = Vim   + (tIM >= cfg.tauIM);
end

need = ceil(cfg.alpha * cfg.Lsym);

out.accept_perm = (Vperm >= need);
out.accept_im   = (Vim   >= need);

end

% ==========================================================
function out = run_one_trial_replay_dppa(cfg, nodes, u0, SNRref_dB, nonce_old, nonce_new)

L = cfg.Lsym;
B = cfg.B;
m = cfg.m;
Npilots = cfg.Npilots;

C_seed = randi([0 1], 1, m);

p_base = ones(Npilots,1);
u_meas = zeros(1,L);
phiH1_deg = zeros(1,L);

snr_lin = 10^(SNRref_dB/10);
sigma2 = 1/snr_lin;

for ell = 1:L
    % TX uses OLD nonce
    nonce_bits_tx = nonce_with_symbol(nonce_old, ell, 16);
    C_eff_tx      = mix_challenge_nonce(C_seed, nonce_bits_tx, m);
    R_tx          = get_puf_bits(nodes(u0), C_eff_tx, B);
    seed_tx       = hash_bits_u32(R_tx);
    phi_tx_deg    = phase_from_seed_dppa(seed_tx, cfg.phi_choices_deg);
    phi_tx        = deg2rad(phi_tx_deg);

    p_tx = p_base * exp(1j*phi_tx);

    w = sqrt(sigma2/2) * (randn(Npilots,1) + 1j*randn(Npilots,1));
    z = p_tx + w;

    u_meas(ell) = (p_base' * z) / (norm(p_base)^2 + 1e-12);

    % Verifier expects NEW nonce
    nonce_bits_H1 = nonce_with_symbol(nonce_new, ell, 16);
    C_eff_H1      = mix_challenge_nonce(C_seed, nonce_bits_H1, m);
    R_H1          = get_puf_bits(nodes(u0), C_eff_H1, B);
    seed_H1       = hash_bits_u32(R_H1);
    phiH1_deg(ell)= phase_from_seed_dppa(seed_H1, cfg.phi_choices_deg);
end

Td_H1 = zeros(1, L-1);

for ell = 2:L
    d_meas = u_meas(ell) * conj(u_meas(ell-1));
    dphi_meas = angle(d_meas);

    dphi_H1 = deg2rad(phiH1_deg(ell) - phiH1_deg(ell-1));
    e1 = angle(exp(1j*(dphi_meas - dphi_H1)));

    Td_H1(ell-1) = cos(e1);
end

Ldiff = numel(Td_H1);
Vd_H1 = sum(Td_H1 >= cfg.tauD);

out.accept_dppa = (Vd_H1 >= ceil(cfg.alpha * Ldiff));

end