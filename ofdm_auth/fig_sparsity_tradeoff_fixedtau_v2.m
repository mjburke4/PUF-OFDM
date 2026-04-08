function fig_sparsity_tradeoff_fixedtau_v2()

cfg = cfg_default();
cfg.chan.type = 'awgn';

cfg.Nfft = 64;
cfg.Npilots = 16;          % set this to match the figure/story

% Operating point for this figure
SNRref = 0;                % stress test
cfg.Lsym  = 16;
cfg.alpha = 0.75;

MC = 6000;
MC_cal = 3000;
targetFAR = 0.01;

% Node setup
U = 10; 
m = 64;
nodes = init_node_population(U, m, 'seed', 7);

% Fix distance so SNRref is meaningful
d0 = 10;
for u = 1:U
    nodes(u).distance_m = d0;
end
u0 = 3;

% Sweep active-pilot sparsity
K_list = 1:(cfg.Npilots-1);
Hbits = zeros(size(K_list));
for i = 1:numel(K_list)
    Hbits(i) = log2(nchoosek(cfg.Npilots, K_list(i)));
end

% ------------------------------------------------------------
% 1) Calibrate tau ONCE at reference sparsity K_cal
% ------------------------------------------------------------
K_cal = floor(cfg.Npilots/2);   % midpoint / max entropy region
cfg_cal = cfg;
cfg_cal.K_active_pilots = K_cal;

tauP_grid  = linspace(0.00, 0.50, 101);
tauIM_grid = linspace(0.00, 1.00, 121);

% ---------- Calibrate tauP ----------
best_tauP = tauP_grid(end);
best_FRR_P = inf;
best_FAR_P = inf;

for tau = tauP_grid
    cfgTmp = cfg_cal;
    cfgTmp.tauP = tau;

    % FAR on H0
    acc = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNRref, false, 0);
        acc = acc + o.accept_perm;
    end
    FAR_here = acc / MC_cal;

    % FRR on H1
    rej = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNRref, true, 0);
        rej = rej + (~o.accept_perm);
    end
    FRR_here = rej / MC_cal;

    if (FAR_here <= targetFAR) && (FRR_here < best_FRR_P)
        best_tauP = tau;
        best_FRR_P = FRR_here;
        best_FAR_P = FAR_here;
    end
end

if isinf(best_FRR_P)
    warning('No tauP met FAR target; using approximate EER-like fallback.');
    best_gap = inf;
    for tau = tauP_grid
        cfgTmp = cfg_cal;
        cfgTmp.tauP = tau;

        acc = 0; rej = 0;
        for t = 1:MC_cal
            o0 = run_one_trial_auth(cfgTmp, nodes, u0, SNRref, false, 0);
            o1 = run_one_trial_auth(cfgTmp, nodes, u0, SNRref, true, 0);
            acc = acc + o0.accept_perm;
            rej = rej + (~o1.accept_perm);
        end
        FAR_here = acc / MC_cal;
        FRR_here = rej / MC_cal;
        gap = abs(FAR_here - FRR_here);

        if gap < best_gap
            best_gap = gap;
            best_tauP = tau;
            best_FRR_P = FRR_here;
            best_FAR_P = FAR_here;
        end
    end
end

% ---------- Calibrate tauIM ----------
best_tauIM = tauIM_grid(end);
best_FRR_IM = inf;
best_FAR_IM = inf;

for tau = tauIM_grid
    cfgTmp = cfg_cal;
    cfgTmp.tauIM = tau;

    % FAR on H0
    acc = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNRref, false, 0);
        acc = acc + o.accept_im;
    end
    FAR_here = acc / MC_cal;

    % FRR on H1
    rej = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNRref, true, 0);
        rej = rej + (~o.accept_im);
    end
    FRR_here = rej / MC_cal;

    if (FAR_here <= targetFAR) && (FRR_here < best_FRR_IM)
        best_tauIM = tau;
        best_FRR_IM = FRR_here;
        best_FAR_IM = FAR_here;
    end
end

if isinf(best_FRR_IM)
    warning('No tauIM met FAR target; using approximate EER-like fallback.');
    best_gap = inf;
    for tau = tauIM_grid
        cfgTmp = cfg_cal;
        cfgTmp.tauIM = tau;

        acc = 0; rej = 0;
        for t = 1:MC_cal
            o0 = run_one_trial_auth(cfgTmp, nodes, u0, SNRref, false, 0);
            o1 = run_one_trial_auth(cfgTmp, nodes, u0, SNRref, true, 0);
            acc = acc + o0.accept_im;
            rej = rej + (~o1.accept_im);
        end
        FAR_here = acc / MC_cal;
        FRR_here = rej / MC_cal;
        gap = abs(FAR_here - FRR_here);

        if gap < best_gap
            best_gap = gap;
            best_tauIM = tau;
            best_FRR_IM = FRR_here;
            best_FAR_IM = FAR_here;
        end
    end
end

cfg.tauP  = best_tauP;
cfg.tauIM = best_tauIM;

fprintf('Global tau calibration at Np=%d, K_cal=%d (H=%.2f bits), SNR=%g dB\n', ...
    cfg.Npilots, K_cal, log2(nchoosek(cfg.Npilots,K_cal)), SNRref);
fprintf('  tauP=%.4f (FAR=%.4f, FRR=%.4f)\n', cfg.tauP, best_FAR_P, best_FRR_P);
fprintf('  tauIM=%.4f (FAR=%.4f, FRR=%.4f)\n\n', cfg.tauIM, best_FAR_IM, best_FRR_IM);

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
    accP = 0; 
    accI = 0;
    for t = 1:MC
        o = run_one_trial_auth(cfgK, nodes, u0, SNRref, false, 0);
        accP = accP + o.accept_perm;
        accI = accI + o.accept_im;
    end
    FARp(ki) = accP / MC;
    FARi(ki) = accI / MC;

    % FRR (H1)
    rejP = 0; 
    rejI = 0;
    for t = 1:MC
        o = run_one_trial_auth(cfgK, nodes, u0, SNRref, true, 0);
        rejP = rejP + (~o.accept_perm);
        rejI = rejI + (~o.accept_im);
    end
    FRRp(ki) = rejP / MC;
    FRRi(ki) = rejI / MC;

    fprintf('K=%2d (H=%.2f bits): FAR(P)=%.4g FRR(P)=%.4g | FAR(IM)=%.4g FRR(IM)=%.4g\n', ...
        K, Hbits(ki), FARp(ki), FRRp(ki), FARi(ki), FRRi(ki));
end

% ------------------------------------------------------------
% 3) Plot FAR/FRR vs K
% ------------------------------------------------------------
figure('Color','w'); 
hold on; grid on; box on;
set(gca,'YScale','log', 'FontSize',12, 'LineWidth',1.0, ...
    'XMinorGrid','off', 'YMinorGrid','on');

plot(K_list, max(FARp,1/MC), '-o', 'LineWidth', 2, 'MarkerSize', 7);
plot(K_list, max(FRRp,1/MC), '--o', 'LineWidth', 2, 'MarkerSize', 7);
plot(K_list, max(FARi,1/MC), '-s', 'LineWidth', 2, 'MarkerSize', 7);
plot(K_list, max(FRRi,1/MC), '--s', 'LineWidth', 2, 'MarkerSize', 7);

xlabel('Active pilots K');
ylabel('Probability (log scale)');
title(sprintf('Sparsity Tradeoff (Npilots=%d, fixed \\tau from K=%d, SNR=%g dB, L=%d, \\alpha=%.2f)', ...
    cfg.Npilots, K_cal, SNRref, cfg.Lsym, cfg.alpha));
legend('FAR Perm','FRR Perm','FAR IM','FRR IM','Location','best');

xlim([min(K_list) max(K_list)]);
ylim([1/MC 1]);

% ------------------------------------------------------------
% 4) Optional entropy plot
% ------------------------------------------------------------
figure('Color','w'); 
grid on; box on;
set(gca,'FontSize',12, 'LineWidth',1.0);

plot(K_list, Hbits, '-d', 'LineWidth', 2, 'MarkerSize', 7);
xlabel('Active pilots K');
ylabel('Mask entropy per symbol: log_2 C(N_p, K) [bits]');
title(sprintf('Combinatorial Mask Entropy (Npilots=%d)', cfg.Npilots));

end