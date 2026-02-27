function fig_far_vs_L()

cfg = cfg_default();
cfg.chan.type = 'awgn';

% Sweep settings
SNR_list = [-5 5 15];
L_list   = [4 8 16 32 64];      % you can trim if runtime gets annoying
MC       = 3000;

% Vote rule: keep alpha fixed for this figure
cfg.alpha = 0.75;

% Calibration settings (single calibration at SNR_cal)
SNR_cal   = 15;
MC_cal    = 2000;
targetFAR = 0.01;              % design target at SNR_cal
Vmax_frac = 0.40;              % constraint on mean H0 votes, fraction of L

% Node setup
U = 10; m = 64;
nodes = init_node_population(U, m, 'seed', 7);

% Fix distance so SNRref is true received SNR at d0
d0 = 10;
for u = 1:U
    nodes(u).distance_m = d0;
end
u0 = 3;

% Preallocate results: [snr_idx x L_idx]
FARp = zeros(numel(SNR_list), numel(L_list));
FRRp = zeros(numel(SNR_list), numel(L_list));
FARi = zeros(numel(SNR_list), numel(L_list));
FRRi = zeros(numel(SNR_list), numel(L_list));

% We'll store calibrated taus per L (since Vmax scales with L)
tauP_cal = zeros(1, numel(L_list));
tauIM_cal = zeros(1, numel(L_list));

for li = 1:numel(L_list)
    cfgL = cfg;
    cfgL.Lsym = L_list(li);

    need = ceil(cfgL.alpha * cfgL.Lsym);
    Vmax = floor(Vmax_frac * cfgL.Lsym);

    % Tau grids (nonnegative)
    tauP_grid  = linspace(0.00, 0.50, 50);
    tauIM_grid = linspace(0.00, 0.50, 50);

    % --- Calibrate tauP at SNR_cal with FAR target + mean-vote constraint
    best_tauP = tauP_grid(end);
    bestFARp = NaN; bestVmeanP = NaN;

    for tau = tauP_grid
        cfgTmp = cfgL; cfgTmp.tauP = tau;

        acc = 0; Vsum = 0;
        for t = 1:MC_cal
            o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, false);
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

    % --- Calibrate tauIM at SNR_cal with FAR target + mean-vote constraint
    best_tauIM = tauIM_grid(end);
    bestFARi = NaN; bestVmeanI = NaN;

    for tau = tauIM_grid
        cfgTmp = cfgL; cfgTmp.tauIM = tau;

        acc = 0; Vsum = 0;
        for t = 1:MC_cal
            o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, false);
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

    cfgL.tauP  = best_tauP;
    cfgL.tauIM = best_tauIM;

    tauP_cal(li)  = best_tauP;
    tauIM_cal(li) = best_tauIM;

    fprintf('L=%d (need %d/%d), calibrated @%gdB: tauP=%.3f (FAR~%.3f, Vmean=%.2f), tauIM=%.3f (FAR~%.3f, Vmean=%.2f)\n', ...
        cfgL.Lsym, need, cfgL.Lsym, SNR_cal, best_tauP, bestFARp, bestVmeanP, best_tauIM, bestFARi, bestVmeanI);

    % --- Now evaluate FAR/FRR vs SNR for this L using calibrated taus
    for si = 1:numel(SNR_list)
        snr = SNR_list(si);

        % Perm: FAR (H0)
        acc = 0;
        for t = 1:MC
            o = run_one_trial_auth(cfgL, nodes, u0, snr, false);
            acc = acc + o.accept_perm;
        end
        FARp(si, li) = acc/MC;

        % Perm: FRR (H1)
        rej = 0;
        for t = 1:MC
            o = run_one_trial_auth(cfgL, nodes, u0, snr, true);
            rej = rej + (~o.accept_perm);
        end
        FRRp(si, li) = rej/MC;

        % IM: FAR (H0)
        acc = 0;
        for t = 1:MC
            o = run_one_trial_auth(cfgL, nodes, u0, snr, false);
            acc = acc + o.accept_im;
        end
        FARi(si, li) = acc/MC;

        % IM: FRR (H1)
        rej = 0;
        for t = 1:MC
            o = run_one_trial_auth(cfgL, nodes, u0, snr, true);
            rej = rej + (~o.accept_im);
        end
        FRRi(si, li) = rej/MC;

        fprintf('  SNR=%g dB: FAR(P)=%.4f FRR(P)=%.4f | FAR(IM)=%.4f FRR(IM)=%.4f\n', ...
            snr, FARp(si,li), FRRp(si,li), FARi(si,li), FRRi(si,li));
    end

end

%% ---------- Plot Permutation ----------
figure; hold on; grid on; set(gca,'YScale','log');
for si = 1:numel(SNR_list)
    plot(L_list, max(FARp(si,:), 1/MC), '-o', 'LineWidth', 3);
    plot(L_list, max(FRRp(si,:), 1/MC), '--s', 'LineWidth', 3);
end
xlabel('Frame length L');
ylabel('Probability (log scale)');
title(sprintf('Permutation: FAR/FRR vs L (\\alpha=%.2f, MC=%d, cal@%gdB)', cfg.alpha, MC, SNR_cal));

leg = strings(1, 2*numel(SNR_list));
for si = 1:numel(SNR_list)
    leg(2*si-1) = sprintf('FAR (SNR=%g dB)', SNR_list(si));
    leg(2*si)   = sprintf('FRR (SNR=%g dB)', SNR_list(si));
end
legend(leg, 'Location','best');

%% ---------- Plot IM ----------
figure; hold on; grid on; set(gca,'YScale','log');
for si = 1:numel(SNR_list)
    plot(L_list, max(FARi(si,:), 1/MC), '-o', 'LineWidth', 3);
    plot(L_list, max(FRRi(si,:), 1/MC), '--s', 'LineWidth', 3);
end
xlabel('Frame length L');
ylabel('Probability (log scale)');
title(sprintf('IM: FAR/FRR vs L (\\alpha=%.2f, MC=%d, cal@%gdB)', cfg.alpha, MC, SNR_cal));

leg = strings(1, 2*numel(SNR_list));
for si = 1:numel(SNR_list)
    leg(2*si-1) = sprintf('FAR (SNR=%g dB)', SNR_list(si));
    leg(2*si)   = sprintf('FRR (SNR=%g dB)', SNR_list(si));
end
legend(leg, 'Location','best');

end