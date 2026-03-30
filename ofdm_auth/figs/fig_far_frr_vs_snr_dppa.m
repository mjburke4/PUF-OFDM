function fig_far_frr_vs_snr_dppa()

cfg = cfg_default_dppa();

U = 10;
nodes = init_node_population(U, cfg.m, 'seed', 7);

u0 = 3;
SNRvec = -5:3:25;
MC = 1000;

% Threshold calibration at one reference SNR, same spirit as perm/IM
SNR_cal = 15;
MC_cal = 400;
targetFAR = 0.01;

tauD_grid = linspace(0.50, 0.999, 120);
best_tauD = cfg.tauD;

Ldiff = cfg.Lsym - 1;
Vmax_frac = 0.40;
Vmax = floor(Vmax_frac * Ldiff);

FAR_est = inf;
Vmean_est = inf;

for tau = tauD_grid
    cfg_tmp = cfg;
    cfg_tmp.tauD = tau;

    acc = 0;
    Vsum = 0;

    for t = 1:MC_cal
        o = run_one_trial_dppa(cfg_tmp, nodes, u0, SNR_cal, false);
        acc  = acc  + o.accept_dppa;
        Vsum = Vsum + o.Vd_H1;   % under impostor TX, verifier still tests H1
    end

    FAR_here = acc / MC_cal;
    Vmean = Vsum / MC_cal;

    if (FAR_here <= targetFAR) && (Vmean <= Vmax)
        best_tauD = tau;
        FAR_est = FAR_here;
        Vmean_est = Vmean;
        break;
    end
end

cfg.tauD = best_tauD;

fprintf('DPPA calibrated at SNR=%g dB:\n', SNR_cal);
fprintf('  tauD=%.4f (FAR~%.3f, mean H0 votes=%.2f)\n', cfg.tauD, FAR_est, Vmean_est);

tauD_list = cfg.tauD * [0.9 0.95 1.00]; %tauD_list = cfg.tauD * [0.95 1 1.05];
tauD_list = min(max(tauD_list, 0), 0.999);

FARd = zeros(numel(tauD_list), numel(SNRvec));
FRRd = zeros(numel(tauD_list), numel(SNRvec));

for si = 1:numel(SNRvec)
    snr = SNRvec(si);

    for ti = 1:numel(tauD_list)
        cfg_tmp = cfg;
        cfg_tmp.tauD = tauD_list(ti);

        % Legit -> FRR
        rej = 0;
        for t = 1:MC
            o = run_one_trial_dppa(cfg_tmp, nodes, u0, snr, true);
            rej = rej + (~o.accept_dppa);
        end
        FRRd(ti, si) = rej / MC;

        % Impostor -> FAR
        acc = 0;
        for t = 1:MC
            o = run_one_trial_dppa(cfg_tmp, nodes, u0, snr, false);
            acc = acc + o.accept_dppa;
        end
        FARd(ti, si) = acc / MC;
    end

    mid = ceil(numel(tauD_list)/2);
    fprintf('SNR=%g dB: FRR(DPPA)=%.3f FAR(DPPA)=%.3f\n', ...
        snr, FRRd(mid,si), FARd(mid,si));
end

figure; hold on; grid on;
for ti = 1:numel(tauD_list)
    plot(SNRvec, FARd(ti,:), '-o', 'LineWidth', 2.5);
    plot(SNRvec, FRRd(ti,:), '--s', 'LineWidth', 2.5);
end
xlabel('SNR_{ref} (dB)');
ylabel('Probability');
title('DPPA Authentication: FAR/FRR vs SNR');

leg = strings(1, 2*numel(tauD_list));
for ti = 1:numel(tauD_list)
    leg(2*ti-1) = sprintf('FAR (\\tau_D=%.3f)', tauD_list(ti));
    leg(2*ti)   = sprintf('FRR (\\tau_D=%.3f)', tauD_list(ti));
end
legend(leg, 'Location','best');

end