function fig_far_frr_vs_L_dppa()

cfg = cfg_default_dppa();
cfg.chan.type = 'awgn';

% Sweep settings
SNR_list = [-5 5 15];
L_list   = [4 8 16 32 64];
MC       = 5000;

% Vote rule
cfg.alpha = 0.75;

% Single calibration settings
L_cal     = 8;      % calibrate tau once at this L
SNR_cal   = 15;
MC_cal    = 3000;
targetFAR = 0.01;
Vmax_frac = 0.40;   % mean H0 vote constraint as fraction of (Ldiff = L-1)

% Node setup
U = 10;
nodes = init_node_population(U, cfg.m, 'seed', 7);
u0 = 3;

% Preallocate
FARd = zeros(numel(SNR_list), numel(L_list));
FRRd = zeros(numel(SNR_list), numel(L_list));

%% =========================================================
%  Single calibration of tauD at (SNR_cal, L_cal)
% ==========================================================
cfg_cal = cfg;
cfg_cal.Lsym = L_cal;

Ldiff_cal = cfg_cal.Lsym - 1;
need_cal  = ceil(cfg_cal.alpha * Ldiff_cal);
Vmax      = floor(Vmax_frac * Ldiff_cal);

tauD_grid = linspace(0.60, 1.00, 81);

best_tauD   = tauD_grid(end);
bestFARd    = NaN;
bestVmeanD  = NaN;

for tau = tauD_grid
    cfgTmp = cfg_cal;
    cfgTmp.tauD = tau;

    acc = 0;
    Vsum = 0;

    for t = 1:MC_cal
        o = run_one_trial_dppa(cfgTmp, nodes, u0, SNR_cal, false);
        acc  = acc  + o.accept_dppa;
        Vsum = Vsum + o.Vd_H1;
    end

    FAR_here = acc / MC_cal;
    Vmean    = Vsum / MC_cal;

    if (FAR_here <= targetFAR) && (Vmean <= Vmax)
        best_tauD  = tau;
        bestFARd   = FAR_here;
        bestVmeanD = Vmean;
        break;
    end
end

cfg.tauD = best_tauD;

fprintf('Single DPPA calibration at SNR=%g dB, L=%d (Ldiff=%d, need %d/%d), Vmax=%d:\n', ...
    SNR_cal, L_cal, Ldiff_cal, need_cal, Ldiff_cal, Vmax);
fprintf('  tauD=%.4f (FAR~%.3f, mean H0 votes=%.2f)\n\n', ...
    cfg.tauD, bestFARd, bestVmeanD);

%% =========================================================
%  Sweep L using FIXED tauD (no recalibration)
% ==========================================================
for li = 1:numel(L_list)
    cfgL = cfg;
    cfgL.Lsym = L_list(li);

    Ldiff = cfgL.Lsym - 1;
    need  = ceil(cfgL.alpha * Ldiff);

    fprintf('Evaluating DPPA L=%d (Ldiff=%d, need %d/%d) with fixed tauD=%.3f\n', ...
        cfgL.Lsym, Ldiff, need, Ldiff, cfgL.tauD);

    for si = 1:numel(SNR_list)
        snr = SNR_list(si);

        % FAR (H0)
        acc = 0;
        for t = 1:MC
            o = run_one_trial_dppa(cfgL, nodes, u0, snr, false);
            acc = acc + o.accept_dppa;
        end
        FARd(si, li) = acc / MC;

        % FRR (H1)
        rej = 0;
        for t = 1:MC
            o = run_one_trial_dppa(cfgL, nodes, u0, snr, true);
            rej = rej + (~o.accept_dppa);
        end
        FRRd(si, li) = rej / MC;

        fprintf('  SNR=%g dB: FAR(DPPA)=%.4f FRR(DPPA)=%.4f\n', ...
            snr, FARd(si,li), FRRd(si,li));
    end
    fprintf('\n');
end

%% ---------- Plot DPPA ----------
figure; hold on; grid on; set(gca,'YScale','log');

for si = 1:numel(SNR_list)
    plot(L_list, max(FARd(si,:), 1/MC), '-o', 'LineWidth', 2);
    plot(L_list, max(FRRd(si,:), 1/MC), '--s', 'LineWidth', 2);
end

xlabel('Frame length L');
ylabel('Probability (log scale)');
title(sprintf('DPPA: FAR/FRR vs L (\\alpha=%.2f, MC=%d, fixed \\tau_D, cal @ SNR=%gdB,L=%d)', ...
    cfg.alpha, MC, SNR_cal, L_cal));

leg = strings(1, 2*numel(SNR_list));
for si = 1:numel(SNR_list)
    leg(2*si-1) = sprintf('FAR (SNR=%g dB)', SNR_list(si));
    leg(2*si)   = sprintf('FRR (SNR=%g dB)', SNR_list(si));
end
legend(leg, 'Location', 'best');

end