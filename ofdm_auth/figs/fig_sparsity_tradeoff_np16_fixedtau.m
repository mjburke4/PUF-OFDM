function fig_sparsity_tradeoff_np16_fixedtau()

cfg = cfg_default();
cfg.chan.type = 'awgn';

cfg.Nfft = 64;
cfg.Npilots = 8;

% Operating point for this figure
SNRref = 5;          % stress test (try 10 if curves are too noisy)
cfg.Lsym  = 4;
cfg.alpha = 0.9;

MC = 5000;           % use 2000 for fast draft, 5000 for final
MC_cal = 3000;       % calibration trials
targetFAR = 0.01;    % design FAR target at calibration K
Vmax_frac = 0.40;    % mean H0 votes constraint (optional but stabilizes)

% Node setup
U = 10; m = 64;
nodes = init_node_population(U, m, 'seed', 7);

% Fix distance so SNRref is meaningful
d0 = 10;
for u=1:U
    nodes(u).distance_m = d0;
end
u0 = 3;

% Sweep active-pilot sparsity
K_list = 1:(cfg.Npilots-1);
Hbits = zeros(size(K_list));
for i=1:numel(K_list)
    Hbits(i) = log2(nchoosek(cfg.Npilots, K_list(i)));
end

% ------------------------------------------------------------
% 1) Calibrate tau ONCE at reference sparsity K_cal
% ------------------------------------------------------------
K_cal = 4;                         % mid-point maximizes entropy for Npilots=16
cfg_cal = cfg;
cfg_cal.K_active_pilots = K_cal;

need = ceil(cfg_cal.alpha * cfg_cal.Lsym);
Vmax = floor(Vmax_frac * cfg_cal.Lsym);

tauP_grid  = linspace(0.00, 0.50, 101);
tauIM_grid = linspace(0.00, 0.50, 101);

% Calibrate tauP
best_tauP = tauP_grid(end);
bestFARp = NaN; bestVmeanP = NaN;

for tau = tauP_grid
    cfgTmp = cfg_cal; cfgTmp.tauP = tau;
    acc=0; Vsum=0;
    for t=1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNRref, false, 0);
        acc  = acc  + o.accept_perm;
        Vsum = Vsum + o.Vperm;
    end
    FAR_here = acc/MC_cal;
    Vmean = Vsum/MC_cal;

    if (FAR_here <= targetFAR) && (Vmean <= Vmax)
        best_tauP = tau;
        bestFARp = FAR_here;
        bestVmeanP = Vmean;
        break;
    end
end

% Calibrate tauIM
best_tauIM = tauIM_grid(end);
bestFARi = NaN; bestVmeanI = NaN;

for tau = tauIM_grid
    cfgTmp = cfg_cal; cfgTmp.tauIM = tau;
    acc=0; Vsum=0;
    for t=1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNRref, false, 0);
        acc  = acc  + o.accept_im;
        Vsum = Vsum + o.Vim;
    end
    FAR_here = acc/MC_cal;
    Vmean = Vsum/MC_cal;

    if (FAR_here <= targetFAR) && (Vmean <= Vmax)
        best_tauIM = tau;
        bestFARi = FAR_here;
        bestVmeanI = Vmean;
        break;
    end
end

cfg.tauP  = best_tauP;
cfg.tauIM = best_tauIM;

fprintf('Global tau calibration at Np=%d, K_cal=%d (H=%.2f bits), SNR=%g dB\n', ...
    cfg.Npilots, K_cal, log2(nchoosek(cfg.Npilots,K_cal)), SNRref);
fprintf('  tauP=%.4f (FAR~%.3f, mean H0 votes=%.2f)\n', cfg.tauP,  bestFARp, bestVmeanP);
fprintf('  tauIM=%.4f (FAR~%.3f, mean H0 votes=%.2f)\n\n', cfg.tauIM, bestFARi, bestVmeanI);

% ------------------------------------------------------------
% 2) Sweep K with FIXED tauP/tauIM
% ------------------------------------------------------------
FARp = zeros(size(K_list));
FRRp = zeros(size(K_list));
FARi = zeros(size(K_list));
FRRi = zeros(size(K_list));

for ki = 1:numel(K_list)
    K = K_list(ki);
    cfgK = cfg;
    cfgK.K_active_pilots = K;

    % FAR (H0)
    accP=0; accI=0;
    for t=1:MC
        o = run_one_trial_auth(cfgK, nodes, u0, SNRref, false, 0);
        accP = accP + o.accept_perm;
        accI = accI + o.accept_im;
    end
    FARp(ki) = accP/MC;
    FARi(ki) = accI/MC;

    % FRR (H1)
    rejP=0; rejI=0;
    for t=1:MC
        o = run_one_trial_auth(cfgK, nodes, u0, SNRref, true, 0);
        rejP = rejP + (~o.accept_perm);
        rejI = rejI + (~o.accept_im);
    end
    FRRp(ki) = rejP/MC;
    FRRi(ki) = rejI/MC;

    fprintf('K=%2d (H=%.2f bits): FAR(P)=%.4g FRR(P)=%.4g | FAR(IM)=%.4g FRR(IM)=%.4g\n', ...
        K, Hbits(ki), FARp(ki), FRRp(ki), FARi(ki), FRRi(ki));
end

% ------------------------------------------------------------
% 3) Plot FAR/FRR vs K (log scale)
% ------------------------------------------------------------
figure; hold on; grid on; set(gca,'YScale','log');
plot(K_list, max(FARp,1/MC), '-o', 'LineWidth', 2);
plot(K_list, max(FRRp,1/MC), '--o', 'LineWidth', 2);
plot(K_list, max(FARi,1/MC), '-s', 'LineWidth', 2);
plot(K_list, max(FRRi,1/MC), '--s', 'LineWidth', 2);

xlabel('Active pilots K');
ylabel('Probability (log scale)');
title(sprintf('Sparsity Tradeoff (Npilots=%d, fixed \\tau from K=%d, SNR=%g dB, L=%d, \\alpha=%.2f)', ...
    cfg.Npilots, K_cal, SNRref, cfg.Lsym, cfg.alpha));
legend('FAR Perm','FRR Perm','FAR IM','FRR IM','Location','best');

% ------------------------------------------------------------
% 4) Optional: Plot entropy vs K (secondary figure)
% ------------------------------------------------------------
figure; grid on;
plot(K_list, Hbits, '-d', 'LineWidth', 2);
xlabel('Active pilots K');
ylabel('Mask entropy per symbol: log_2 C(N_p, K) [bits]');
title(sprintf('Combinatorial Mask Entropy (Npilots=%d)', cfg.Npilots));

end