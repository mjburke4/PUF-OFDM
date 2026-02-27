function fig_sparsity_tradeoff_np8()

cfg = cfg_default();
cfg.chan.type = 'awgn';

cfg.Nfft = 64;
cfg.Npilots = 8;

SNRref = 5;      % choose 5 or 10 (stress vs moderate)
cfg.Lsym  = 8;
cfg.alpha = 0.75;

MC = 2000;       % bump to 5000 for final
MC_cal = 2000;
targetFAR = 0.01;
Vmax_frac = 0.40;

U = 10; m = 64;
nodes = init_node_population(U, m, 'seed', 7);

d0 = 10;
for u=1:U
    nodes(u).distance_m = d0;
end
u0 = 3;

K_list = 1:(cfg.Npilots-1);

FARp = zeros(size(K_list));
FRRp = zeros(size(K_list));
FARi = zeros(size(K_list));
FRRi = zeros(size(K_list));
Hbits = zeros(size(K_list)); % log2 comb

for ki = 1:numel(K_list)
    K = K_list(ki);
    cfgK = cfg;
    cfgK.K_active_pilots = K;

    % entropy of mask space per symbol (pilot-domain IM)
    Hbits(ki) = log2(nchoosek(cfg.Npilots, K));

    % ---- Calibrate tauP and tauIM ONCE PER K (keeps FAR target consistent)
    % If you want fixed tau across K, tell me and I'll change it.
    need = ceil(cfgK.alpha * cfgK.Lsym);
    Vmax = floor(Vmax_frac * cfgK.Lsym);

    tauP_grid  = linspace(0.00, 0.50, 101);
    tauIM_grid = linspace(0.00, 0.50, 101);

    % Calibrate tauP
    best_tauP = tauP_grid(end);
    for tau = tauP_grid
        cfgTmp = cfgK; cfgTmp.tauP = tau;
        acc=0; Vsum=0;
        for t=1:MC_cal
            o = run_one_trial_auth(cfgTmp, nodes, u0, SNRref, false, 0);
            acc = acc + o.accept_perm;
            Vsum = Vsum + o.Vperm;
        end
        FAR_here = acc/MC_cal;
        Vmean = Vsum/MC_cal;
        if (FAR_here <= targetFAR) && (Vmean <= floor(Vmax_frac*cfgK.Lsym))
            best_tauP = tau; break;
        end
    end

    % Calibrate tauIM
    best_tauIM = tauIM_grid(end);
    for tau = tauIM_grid
        cfgTmp = cfgK; cfgTmp.tauIM = tau;
        acc=0; Vsum=0;
        for t=1:MC_cal
            o = run_one_trial_auth(cfgTmp, nodes, u0, SNRref, false, 0);
            acc = acc + o.accept_im;
            Vsum = Vsum + o.Vim;
        end
        FAR_here = acc/MC_cal;
        Vmean = Vsum/MC_cal;
        if (FAR_here <= targetFAR) && (Vmean <= floor(Vmax_frac*cfgK.Lsym))
            best_tauIM = tau; break;
        end
    end

    cfgK.tauP  = best_tauP;
    cfgK.tauIM = best_tauIM;

    % ---- Evaluate FAR/FRR at this K
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

    fprintf('Np=%d K=%d (H=%.2f bits): FAR(P)=%.3g FRR(P)=%.3g | FAR(IM)=%.3g FRR(IM)=%.3g\n', ...
        cfg.Npilots, K, Hbits(ki), FARp(ki), FRRp(ki), FARi(ki), FRRi(ki));
end

% ---- Plot
figure; hold on; grid on; set(gca,'YScale','log');
plot(K_list, max(FARp,1/MC), '-o', 'LineWidth', 2);
plot(K_list, max(FRRp,1/MC), '--o', 'LineWidth', 2);
plot(K_list, max(FARi,1/MC), '-s', 'LineWidth', 2);
plot(K_list, max(FRRi,1/MC), '--s', 'LineWidth', 2);

xlabel('Active pilots K');
ylabel('Probability (log scale)');
title(sprintf('Sparsity Tradeoff (Npilots=%d, SNR=%g dB, L=%d, \\alpha=%.2f)', ...
    cfg.Npilots, SNRref, cfg.Lsym, cfg.alpha));

legend('FAR Perm','FRR Perm','FAR IM','FRR IM','Location','best');

end