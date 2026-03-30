function fig_eer_vs_snr_compare_perm_im_dppa()

%% ---------------- Common setup ----------------
cfg = cfg_default();                 % your perm/IM cfg
cfg.chan.type = 'awgn';
cfg.Lsym  = 8;
cfg.alpha = 0.75;

cfgD = cfg_default_dppa();          % your DPPA cfg
cfgD.Lsym  = cfg.Lsym;
cfgD.alpha = cfg.alpha;

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
MC = 1500;   % can start with 500 if slow, then increase

% Threshold grids
tauP_grid  = linspace(0, 1, 121);
tauIM_grid = linspace(0, 1, 121);
tauD_grid  = linspace(0.5, 0.999, 121);

EER_perm = zeros(size(SNRvec));
EER_im   = zeros(size(SNRvec));
EER_dppa = zeros(size(SNRvec));

tauP_star  = zeros(size(SNRvec));
tauIM_star = zeros(size(SNRvec));
tauD_star  = zeros(size(SNRvec));

%% ---------------- Main SNR sweep ----------------
for si = 1:numel(SNRvec)
    snr = SNRvec(si);
    fprintf('\nSNR = %g dB\n', snr);

    %% ===== Permutation EER =====
    FARp = zeros(size(tauP_grid));
    FRRp = zeros(size(tauP_grid));

    for ti = 1:numel(tauP_grid)
        cfg.tauP = tauP_grid(ti);

        rej = 0;
        acc = 0;
        for t = 1:MC
            o1 = run_one_trial_auth(cfg, nodes, u0, snr, true);   % legit
            o0 = run_one_trial_auth(cfg, nodes, u0, snr, false);  % impostor
            rej = rej + (~o1.accept_perm);
            acc = acc + (o0.accept_perm);
        end
        FRRp(ti) = rej / MC;
        FARp(ti) = acc / MC;
    end

    [EER_perm(si), tauP_star(si)] = estimate_eer_from_curves(tauP_grid, FARp, FRRp);

    %% ===== IM EER =====
    FARi = zeros(size(tauIM_grid));
    FRRi = zeros(size(tauIM_grid));

    for ti = 1:numel(tauIM_grid)
        cfg.tauIM = tauIM_grid(ti);

        rej = 0;
        acc = 0;
        for t = 1:MC
            o1 = run_one_trial_auth(cfg, nodes, u0, snr, true);
            o0 = run_one_trial_auth(cfg, nodes, u0, snr, false);
            rej = rej + (~o1.accept_im);
            acc = acc + (o0.accept_im);
        end
        FRRi(ti) = rej / MC;
        FARi(ti) = acc / MC;
    end

    [EER_im(si), tauIM_star(si)] = estimate_eer_from_curves(tauIM_grid, FARi, FRRi);

    %% ===== DPPA EER =====
    FARd = zeros(size(tauD_grid));
    FRRd = zeros(size(tauD_grid));

    for ti = 1:numel(tauD_grid)
        cfgD.tauD = tauD_grid(ti);

        rej = 0;
        acc = 0;
        for t = 1:MC
            o1 = run_one_trial_dppa(cfgD, nodes, u0, snr, true);
            o0 = run_one_trial_dppa(cfgD, nodes, u0, snr, false);

            rej = rej + (~o1.accept_dppa);
            acc = acc + (o0.accept_dppa);
        end
        FRRd(ti) = rej / MC;
        FARd(ti) = acc / MC;
    end

    [EER_dppa(si), tauD_star(si)] = estimate_eer_from_curves(tauD_grid, FARd, FRRd);

    fprintf('  Perm: EER=%.4f at tauP=%.3f\n',  EER_perm(si), tauP_star(si));
    fprintf('  IM  : EER=%.4f at tauIM=%.3f\n', EER_im(si),   tauIM_star(si));
    fprintf('  DPPA: EER=%.4f at tauD=%.3f\n',  EER_dppa(si), tauD_star(si));
end
plot_floor = 1e-5;
EER_perm_plot = max(EER_perm, plot_floor);
EER_im_plot   = max(EER_im,   plot_floor);
EER_dppa_plot = max(EER_dppa, plot_floor);

%% ---------------- Comparison plot ----------------
figure; hold on; grid on;
plot(SNRvec, EER_perm, '-o', 'LineWidth', 2);
plot(SNRvec, EER_im,   '-s', 'LineWidth', 2);
plot(SNRvec, EER_dppa, '-^', 'LineWidth', 2);

xlabel('SNR_{ref} at d_0 (dB)');
ylabel('Equal Error Rate (EER)');
title(sprintf('Authentication Comparison: EER vs SNR (L=%d, \\alpha=%.2f)', ...
    cfg.Lsym, cfg.alpha));
set(gca, 'YScale', 'log');
legend('Permutation', 'IM', 'DPPA', 'Location', 'best');

%% ---------------- Optional: threshold trajectories ----------------
figure; hold on; grid on;
plot(SNRvec, tauP_star,  '-o', 'LineWidth', 2);
plot(SNRvec, tauIM_star, '-s', 'LineWidth', 2);
plot(SNRvec, tauD_star,  '-^', 'LineWidth', 2);
xlabel('SNR_{ref} at d_0 (dB)');
ylabel('Threshold at EER');
title('Threshold yielding EER vs SNR');
legend('\tau_P^*', '\tau_{IM}^*', '\tau_D^*', 'Location', 'best');

end

function [eer, tau_star] = estimate_eer_from_curves(tau_grid, FAR, FRR)
% Estimate EER by choosing the threshold where |FAR - FRR| is minimized

    [~, idx] = min(abs(FAR - FRR));
    eer = 0.5 * (FAR(idx) + FRR(idx));
    tau_star = tau_grid(idx);
end