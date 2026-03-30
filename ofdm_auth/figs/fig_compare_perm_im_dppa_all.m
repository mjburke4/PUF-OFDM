function fig_compare_perm_im_dppa_all()

%% ---------------- Common setup ----------------
cfg = cfg_default();                 % Perm/IM cfg
cfg.chan.type = 'awgn';
%cfg.chan.type = 'rayleigh';
cfg.Lsym  = 16;
cfg.alpha = 0.75;

tic

cfgD = cfg_default_dppa();          % DPPA cfg
cfgD.Lsym  = cfg.Lsym;
cfgD.alpha = cfg.alpha;

cfg.adv_mode  = 'impostor_othernode';%'; %impostor_othernode';   % or 'impostor_random'
cfgD.adv_mode = cfg.adv_mode;

U = 10;
m = 64;
nodes = init_node_population(U, m, 'seed', 7);

% Fix all nodes at reference distance so SNRref is true SNR
d0 = 10;
for u = 1:U
    nodes(u).distance_m = d0;
end

u0 = 3;
SNRvec = -5:3:25;

% Monte Carlo controls
MC_eer = 1500;       % for EER sweeps
MC_cal = 2000;       % for calibration at reference SNR
MC_op  = 1500;       % for fixed-threshold FAR/FRR vs SNR

% Calibration settings
SNR_cal = 15;
targetFAR = 1e-2;

% Threshold grids
tauP_grid  = linspace(0, 1, 100);
tauIM_grid = linspace(0, 1, 100);
tauD_grid  = linspace(0, 0.999, 100); %linspace(0.5, 0.999, 100);

%% ---------------- Storage ----------------
EER_perm = zeros(size(SNRvec));
EER_im   = zeros(size(SNRvec));
EER_dppa = zeros(size(SNRvec));

tauP_star  = zeros(size(SNRvec));
tauIM_star = zeros(size(SNRvec));
tauD_star  = zeros(size(SNRvec));

% Fixed-threshold operating curves
FAR_perm_fix = zeros(size(SNRvec));
FRR_perm_fix = zeros(size(SNRvec));

FAR_im_fix   = zeros(size(SNRvec));
FRR_im_fix   = zeros(size(SNRvec));

FAR_dppa_fix = zeros(size(SNRvec));
FRR_dppa_fix = zeros(size(SNRvec));

%% ==========================================================
% 1) Threshold calibration at SNR_cal for target FAR
%% ==========================================================
fprintf('\n========================================\n');
fprintf('Threshold calibration at SNR = %g dB\n', SNR_cal);
fprintf('Target FAR <= %.3g\n', targetFAR);
fprintf('========================================\n');

% ---- Perm calibration
tauP_cal = tauP_grid(end);
for ti = 1:numel(tauP_grid)
    cfg.tauP = tauP_grid(ti);

    acc = 0;
    for t = 1:MC_cal
        o0 = run_one_trial_auth(cfg, nodes, u0, SNR_cal, false);
        acc = acc + o0.accept_perm;
    end
    FAR_here = acc / MC_cal;

    if FAR_here <= targetFAR
        tauP_cal = tauP_grid(ti);
        break;
    end
end

% ---- IM calibration
tauIM_cal = tauIM_grid(end);
for ti = 1:numel(tauIM_grid)
    cfg.tauIM = tauIM_grid(ti);

    acc = 0;
    for t = 1:MC_cal
        o0 = run_one_trial_auth(cfg, nodes, u0, SNR_cal, false);
        acc = acc + o0.accept_im;
    end
    FAR_here = acc / MC_cal;

    if FAR_here <= targetFAR
        tauIM_cal = tauIM_grid(ti);
        break;
    end
end

% ---- DPPA calibration
tauD_cal = tauD_grid(end);
for ti = 1:numel(tauD_grid)
    cfgD.tauD = tauD_grid(ti);

    acc = 0;
    for t = 1:MC_cal
        o0 = run_one_trial_dppa(cfgD, nodes, u0, SNR_cal, false);
        acc = acc + o0.accept_dppa;
    end
    FAR_here = acc / MC_cal;

    if FAR_here <= targetFAR
        tauD_cal = tauD_grid(ti);
        break;
    end
end

fprintf('Calibrated thresholds at %g dB:\n', SNR_cal);
fprintf('  tauP_cal  = %.4f\n', tauP_cal);
fprintf('  tauIM_cal = %.4f\n', tauIM_cal);
fprintf('  tauD_cal  = %.4f\n', tauD_cal);

%% ==========================================================
% 2) Main SNR sweep: EER + fixed-threshold operational curves
%% ==========================================================
for si = 1:numel(SNRvec)
    snr = SNRvec(si);
    fprintf('\nSNR = %g dB\n', snr);

    %% ===== Perm EER =====
    FARp = zeros(size(tauP_grid));
    FRRp = zeros(size(tauP_grid));

    for ti = 1:numel(tauP_grid)
        cfg.tauP = tauP_grid(ti);

        rej = 0;
        acc = 0;
        for t = 1:MC_eer
            o1 = run_one_trial_auth(cfg, nodes, u0, snr, true);
            o0 = run_one_trial_auth(cfg, nodes, u0, snr, false);

            rej = rej + (~o1.accept_perm);
            acc = acc + (o0.accept_perm);
        end
        FRRp(ti) = rej / MC_eer;
        FARp(ti) = acc / MC_eer;
    end

    [EER_perm(si), tauP_star(si)] = estimate_eer_from_curves(tauP_grid, FARp, FRRp);

    %% ===== IM EER =====
    FARi = zeros(size(tauIM_grid));
    FRRi = zeros(size(tauIM_grid));

    for ti = 1:numel(tauIM_grid)
        cfg.tauIM = tauIM_grid(ti);

        rej = 0;
        acc = 0;
        for t = 1:MC_eer
            o1 = run_one_trial_auth(cfg, nodes, u0, snr, true);
            o0 = run_one_trial_auth(cfg, nodes, u0, snr, false);

            rej = rej + (~o1.accept_im);
            acc = acc + (o0.accept_im);
        end
        FRRi(ti) = rej / MC_eer;
        FARi(ti) = acc / MC_eer;
    end

    [EER_im(si), tauIM_star(si)] = estimate_eer_from_curves(tauIM_grid, FARi, FRRi);

    %% ===== DPPA EER =====
    FARd = zeros(size(tauD_grid));
    FRRd = zeros(size(tauD_grid));

    for ti = 1:numel(tauD_grid)
        cfgD.tauD = tauD_grid(ti);

        rej = 0;
        acc = 0;
        for t = 1:MC_eer
            o1 = run_one_trial_dppa(cfgD, nodes, u0, snr, true);
            o0 = run_one_trial_dppa(cfgD, nodes, u0, snr, false);

            rej = rej + (~o1.accept_dppa);
            acc = acc + (o0.accept_dppa);
        end
        FRRd(ti) = rej / MC_eer;
        FARd(ti) = acc / MC_eer;
    end

    [EER_dppa(si), tauD_star(si)] = estimate_eer_from_curves(tauD_grid, FARd, FRRd);

    %% ===== Perm fixed-threshold operational performance =====
    cfg.tauP = tauP_cal;
    rej = 0; acc = 0;
    for t = 1:MC_op
        o1 = run_one_trial_auth(cfg, nodes, u0, snr, true);
        o0 = run_one_trial_auth(cfg, nodes, u0, snr, false);
        rej = rej + (~o1.accept_perm);
        acc = acc + (o0.accept_perm);
    end
    FRR_perm_fix(si) = rej / MC_op;
    FAR_perm_fix(si) = acc / MC_op;

    %% ===== IM fixed-threshold operational performance =====
    cfg.tauIM = tauIM_cal;
    rej = 0; acc = 0;
    for t = 1:MC_op
        o1 = run_one_trial_auth(cfg, nodes, u0, snr, true);
        o0 = run_one_trial_auth(cfg, nodes, u0, snr, false);
        rej = rej + (~o1.accept_im);
        acc = acc + (o0.accept_im);
    end
    FRR_im_fix(si) = rej / MC_op;
    FAR_im_fix(si) = acc / MC_op;

    %% ===== DPPA fixed-threshold operational performance =====
    cfgD.tauD = tauD_cal;
    rej = 0; acc = 0;
    for t = 1:MC_op
        o1 = run_one_trial_dppa(cfgD, nodes, u0, snr, true);
        o0 = run_one_trial_dppa(cfgD, nodes, u0, snr, false);
        rej = rej + (~o1.accept_dppa);
        acc = acc + (o0.accept_dppa);
    end
    FRR_dppa_fix(si) = rej / MC_op;
    FAR_dppa_fix(si) = acc / MC_op;

    fprintf('  Perm: EER=%.4f at tauP=%.3f | FRR=%.4f FAR=%.4f at tauP_cal=%.3f\n', ...
        EER_perm(si), tauP_star(si), FRR_perm_fix(si), FAR_perm_fix(si), tauP_cal);
    fprintf('  IM  : EER=%.4f at tauIM=%.3f | FRR=%.4f FAR=%.4f at tauIM_cal=%.3f\n', ...
        EER_im(si), tauIM_star(si), FRR_im_fix(si), FAR_im_fix(si), tauIM_cal);
    fprintf('  DPPA: EER=%.4f at tauD=%.3f | FRR=%.4f FAR=%.4f at tauD_cal=%.3f\n', ...
        EER_dppa(si), tauD_star(si), FRR_dppa_fix(si), FAR_dppa_fix(si), tauD_cal);
end

%% ---------------- Floors for log plotting ----------------
plot_floor = 1e-5;

EER_perm_plot = max(EER_perm, plot_floor);
EER_im_plot   = max(EER_im,   plot_floor);
EER_dppa_plot = max(EER_dppa, plot_floor);

FAR_perm_fix_plot = max(FAR_perm_fix, plot_floor);
FRR_perm_fix_plot = max(FRR_perm_fix, plot_floor);

FAR_im_fix_plot   = max(FAR_im_fix, plot_floor);
FRR_im_fix_plot   = max(FRR_im_fix, plot_floor);

FAR_dppa_fix_plot = max(FAR_dppa_fix, plot_floor);
FRR_dppa_fix_plot = max(FRR_dppa_fix, plot_floor);

%% ==========================================================
% 3) Plot: EER vs SNR
%% ==========================================================
figure; hold on; grid on;
plot(SNRvec, EER_perm_plot, '-o', 'LineWidth', 2);
plot(SNRvec, EER_im_plot,   '-s', 'LineWidth', 2);
plot(SNRvec, EER_dppa_plot, '-^', 'LineWidth', 2);

xlabel('SNR_{ref} at d_0 (dB)');
ylabel('Equal Error Rate (EER)');
title(sprintf('Authentication Comparison: EER vs SNR (L=%d, \\alpha=%.2f)', ...
    cfg.Lsym, cfg.alpha));
set(gca, 'YScale', 'log');
legend('Permutation', 'IM', 'DPPA', 'Location', 'best');

%% ==========================================================
% 4) Plot: threshold yielding EER vs SNR
%% ==========================================================
figure; hold on; grid on;
plot(SNRvec, tauP_star,  '-o', 'LineWidth', 2);
plot(SNRvec, tauIM_star, '-s', 'LineWidth', 2);
plot(SNRvec, tauD_star,  '-^', 'LineWidth', 2);
xlabel('SNR_{ref} at d_0 (dB)');
ylabel('Threshold at EER');
title('Threshold yielding EER vs SNR');
legend('\tau_P^*', '\tau_{IM}^*', '\tau_D^*', 'Location', 'best');

%% ==========================================================
% 5) Plot: FAR vs SNR with thresholds fixed at 15 dB
%% ==========================================================
figure; hold on; grid on;
plot(SNRvec, FAR_perm_fix_plot, '-o', 'LineWidth', 2);
plot(SNRvec, FAR_im_fix_plot,   '-s', 'LineWidth', 2);
plot(SNRvec, FAR_dppa_fix_plot, '-^', 'LineWidth', 2);

xlabel('SNR_{ref} at d_0 (dB)');
ylabel('False Acceptance Rate (FAR)');
title(sprintf('FAR vs SNR with thresholds calibrated at %d dB', SNR_cal));
set(gca, 'YScale', 'log');
legend( ...
    sprintf('Permutation (\\tau_P=%.3f)', tauP_cal), ...
    sprintf('IM (\\tau_{IM}=%.3f)', tauIM_cal), ...
    sprintf('DPPA (\\tau_D=%.3f)', tauD_cal), ...
    'Location', 'best');

%% ==========================================================
% 6) Plot: FRR vs SNR with thresholds fixed at 15 dB
%% ==========================================================
figure; hold on; grid on;
plot(SNRvec, FRR_perm_fix_plot, '-o', 'LineWidth', 2);
plot(SNRvec, FRR_im_fix_plot,   '-s', 'LineWidth', 2);
plot(SNRvec, FRR_dppa_fix_plot, '-^', 'LineWidth', 2);

xlabel('SNR_{ref} at d_0 (dB)');
ylabel('False Rejection Rate (FRR)');
title(sprintf('FRR vs SNR with thresholds calibrated at %d dB', SNR_cal));
set(gca, 'YScale', 'log');
legend( ...
    sprintf('Permutation (\\tau_P=%.3f)', tauP_cal), ...
    sprintf('IM (\\tau_{IM}=%.3f)', tauIM_cal), ...
    sprintf('DPPA (\\tau_D=%.3f)', tauD_cal), ...
    'Location', 'best');

%% ==========================================================
% 7) Optional combined operational plot
%% ==========================================================
figure; hold on; grid on;
plot(SNRvec, FAR_perm_fix_plot, '-o', 'LineWidth', 2);
plot(SNRvec, FRR_perm_fix_plot, '--o', 'LineWidth', 2);

plot(SNRvec, FAR_im_fix_plot, '-s', 'LineWidth', 2);
plot(SNRvec, FRR_im_fix_plot, '--s', 'LineWidth', 2);

plot(SNRvec, FAR_dppa_fix_plot, '-^', 'LineWidth', 2);
plot(SNRvec, FRR_dppa_fix_plot, '--^', 'LineWidth', 2);

xlabel('SNR_{ref} at d_0 (dB)');
ylabel('Probability');
title(sprintf('Operational FAR/FRR vs SNR (thresholds calibrated at %d dB)', SNR_cal));
set(gca, 'YScale', 'log');
legend( ...
    'Perm FAR', 'Perm FRR', ...
    'IM FAR', 'IM FRR', ...
    'DPPA FAR', 'DPPA FRR', ...
    'Location', 'best');
toc

end

function [eer, tau_star] = estimate_eer_from_curves(tau_grid, FAR, FRR)
% Estimate EER by choosing the threshold where |FAR - FRR| is minimized
    [~, idx] = min(abs(FAR - FRR));
    eer = 0.5 * (FAR(idx) + FRR(idx));
    tau_star = tau_grid(idx);
end