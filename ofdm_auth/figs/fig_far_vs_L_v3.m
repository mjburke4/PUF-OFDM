function fig_far_vs_L()

cfg = cfg_default();
cfg.chan.type = 'awgn';

% Sweep settings
SNR_list = [-5 5 15];
L_list   = [4 8 16 32 64];
MC       = 5000;          % increase if you want smoother tails

% Vote rule
cfg.alpha = 0.75;

% Single calibration settings
L_cal     = 8;            % calibrate tau once at this L
SNR_cal   = 15;
MC_cal    = 3000;
targetFAR = 0.01;         % FAR target at SNR_cal
Vmax_frac = 0.40;         % mean H0 vote constraint (fraction of L_cal)

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

%% =========================================================
%  Single calibration of tauP and tauIM at (SNR_cal, L_cal)
% ==========================================================
cfg_cal = cfg;
cfg_cal.Lsym = L_cal;

need_cal = ceil(cfg_cal.alpha * cfg_cal.Lsym);
Vmax = floor(Vmax_frac * cfg_cal.Lsym);

tauP_grid  = linspace(0.00, 0.50, 101);
tauIM_grid = linspace(0.00, 0.50, 101);

% --- Calibrate tauP once
best_tauP = tauP_grid(end);
bestFARp = NaN; bestVmeanP = NaN;

for tau = tauP_grid
    cfgTmp = cfg_cal; cfgTmp.tauP = tau;

    acc = 0; Vsum = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, false);  % H0
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

% --- Calibrate tauIM once
best_tauIM = tauIM_grid(end);
bestFARi = NaN; bestVmeanI = NaN;

for tau = tauIM_grid
    cfgTmp = cfg_cal; cfgTmp.tauIM = tau;

    acc = 0; Vsum = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, false);  % H0
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

% Freeze calibrated thresholds globally
cfg.tauP  = best_tauP;
cfg.tauIM = best_tauIM;

fprintf('Single calibration at SNR=%g dB, L=%d (need %d/%d), Vmax=%d:\n', ...
    SNR_cal, L_cal, need_cal, L_cal, Vmax);
fprintf('  tauP=%.4f (FAR~%.3f, mean H0 votes=%.2f)\n', cfg.tauP,  bestFARp, bestVmeanP);
fprintf('  tauIM=%.4f (FAR~%.3f, mean H0 votes=%.2f)\n\n', cfg.tauIM, bestFARi, bestVmeanI);

%% =========================================================
%  Sweep L using FIXED tauP/tauIM (no recalibration)
% ==========================================================
for li = 1:numel(L_list)
    cfgL = cfg;
    cfgL.Lsym = L_list(li);

    need = ceil(cfgL.alpha * cfgL.Lsym);

    fprintf('Evaluating L=%d (need %d/%d) with fixed tauP=%.3f, tauIM=%.3f\n', ...
        cfgL.Lsym, need, cfgL.Lsym, cfgL.tauP, cfgL.tauIM);

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

    fprintf('\n');
end

%% ---------- Plot Permutation ----------
figure; hold on; grid on; set(gca,'YScale','log');
for si = 1:numel(SNR_list)
    plot(L_list, max(FARp(si,:), 1/MC), '-o', 'LineWidth', 2);
    plot(L_list, max(FRRp(si,:), 1/MC), '--s', 'LineWidth', 2);
end
xlabel('Frame length L');
ylabel('Probability (log scale)');
title(sprintf('Permutation: FAR/FRR vs L (\\alpha=%.2f, MC=%d, fixed \\tau, cal @ SNR=%gdB,L=%d)', ...
    cfg.alpha, MC, SNR_cal, L_cal));

leg = strings(1, 2*numel(SNR_list));
for si = 1:numel(SNR_list)
    leg(2*si-1) = sprintf('FAR (SNR=%g dB)', SNR_list(si));
    leg(2*si)   = sprintf('FRR (SNR=%g dB)', SNR_list(si));
end
legend(leg, 'Location','best');

%% ---------- Plot IM ----------
figure; hold on; grid on; set(gca,'YScale','log');
for si = 1:numel(SNR_list)
    plot(L_list, max(FARi(si,:), 1/MC), '-o', 'LineWidth', 2);
    plot(L_list, max(FRRi(si,:), 1/MC), '--s', 'LineWidth', 2);
end
xlabel('Frame length L');
ylabel('Probability (log scale)');
title(sprintf('IM: FAR/FRR vs L (\\alpha=%.2f, MC=%d, fixed \\tau, cal @ SNR=%gdB,L=%d)', ...
    cfg.alpha, MC, SNR_cal, L_cal));

leg = strings(1, 2*numel(SNR_list));
for si = 1:numel(SNR_list)
    leg(2*si-1) = sprintf('FAR (SNR=%g dB)', SNR_list(si));
    leg(2*si)   = sprintf('FRR (SNR=%g dB)', SNR_list(si));
end
legend(leg, 'Location','best');

end