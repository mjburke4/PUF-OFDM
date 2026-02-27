% -------- Single tau calibration (once) --------
L_cal   = 8;
SNR_cal = 15;
MC_cal  = 4000;
targetFAR = 0.01;
Vmax_frac = 0.40;

cfg = cfg_default();
cfg_cal = cfg;
cfg_cal.Lsym = L_cal;

need_cal = ceil(cfg_cal.alpha * cfg_cal.Lsym);
Vmax = floor(Vmax_frac * cfg_cal.Lsym);

tauP_grid  = linspace(0.00, 0.50, 101);
tauIM_grid = linspace(0.00, 0.50, 101);

% Calibrate tauP once
best_tauP = tauP_grid(end);
bestFARp = NaN; bestVmeanP = NaN;

for tau = tauP_grid
    cfg_tmp = cfg_cal; cfg_tmp.tauP = tau;
    acc = 0; Vsum = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfg_tmp, nodes, u0, SNR_cal, false);
        acc = acc + o.accept_perm;
        Vsum = Vsum + o.Vperm;
    end
    FAR_here = acc/MC_cal;
    Vmean = Vsum/MC_cal;
    if (FAR_here <= targetFAR) && (Vmean <= Vmax)
        best_tauP = tau; bestFARp = FAR_here; bestVmeanP = Vmean;
        break;
    end
end

% Calibrate tauIM once
best_tauIM = tauIM_grid(end);
bestFARi = NaN; bestVmeanI = NaN;

for tau = tauIM_grid
    cfg_tmp = cfg_cal; cfg_tmp.tauIM = tau;
    acc = 0; Vsum = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfg_tmp, nodes, u0, SNR_cal, false);
        acc = acc + o.accept_im;
        Vsum = Vsum + o.Vim;
    end
    FAR_here = acc/MC_cal;
    Vmean = Vsum/MC_cal;
    if (FAR_here <= targetFAR) && (Vmean <= Vmax)
        best_tauIM = tau; bestFARi = FAR_here; bestVmeanI = Vmean;
        break;
    end
end

cfg.tauP  = best_tauP;
cfg.tauIM = best_tauIM;

fprintf('Single calibration at SNR=%g dB, L=%d (need %d/%d), Vmax=%d:\n', ...
    SNR_cal, L_cal, need_cal, L_cal, Vmax);
fprintf('  tauP=%.4f (FAR~%.3f, mean H0 votes=%.2f)\n', cfg.tauP, bestFARp, bestVmeanP);
fprintf('  tauIM=%.4f (FAR~%.3f, mean H0 votes=%.2f)\n', cfg.tauIM, bestFARi, bestVmeanI);
% ----------------------------------------------