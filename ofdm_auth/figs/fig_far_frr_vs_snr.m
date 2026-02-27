function fig_far_frr_vs_snr()

cfg = cfg_default();
cfg.chan.type = 'awgn';

% IMPORTANT: set thresholds (you can tune later or sweep)
cfg.tauP  = 0.05;   % your current good ballpark at ~-12 dB eff SNR
cfg.tauIM = 0.10;   % set from your observed Tim averages
cfg.alpha = 0.75; %0.9
cfg.Lsym  = 8; %32;

U = 10; m = 64;
nodes = init_node_population(U, m, 'seed', 7);

% --- Fix distance so SNRref means what we think it means
d0 = 10;  % must match run_one_trial_auth
for u = 1:U
    nodes(u).distance_m = d0;
end


u0 = 3;

SNRvec = -5:3:25;        % SNRref at d0
MC = 1000; %300;                % trials per SNR point (start small if slow)

% ---------- Calibrate thresholds at a reference SNR ----------
SNR_cal = 15;      % calibration point (dB)
MC_cal  = 400;     % calibration trials
targetFAR = 0.01;  % aim FAR ~ 1% at SNR_cal

% Collect H0 per-symbol statistics by using votes directly:
Vperm0 = zeros(1, MC_cal);
Vim0   = zeros(1, MC_cal);
Vperm1 = zeros(1, MC_cal);
Vim1   = zeros(1, MC_cal);

for t = 1:MC_cal
    o0 = run_one_trial_auth(cfg, nodes, u0, SNR_cal, false); % impostor
    o1 = run_one_trial_auth(cfg, nodes, u0, SNR_cal, true);  % legit
    Vperm0(t) = o0.Vperm;  Vim0(t) = o0.Vim;
    Vperm1(t) = o1.Vperm;  Vim1(t) = o1.Vim;
end

% We choose per-symbol thresholds indirectly by choosing vote thresholds.
% Since your decision is votes >= ceil(alpha*L), the only knobs left are tauP/tauIM.
% But easiest: set tauP/tauIM so that H0 vote counts are usually below ceil(alpha*L).
% We'll instead tune tauP/tauIM by sweeping candidate thresholds quickly.

tauP_grid  = linspace(0.00, 0.50, 101);   % permutation threshold on real(T)
tauIM_grid = linspace(0.00, 0.50, 101);   % IM normalized threshold (>=0)

best_tauP  = cfg.tauP;
best_tauIM = cfg.tauIM;

L    = cfg.Lsym;
need = ceil(cfg.alpha * L);

% Secondary constraint: keep H0 mean votes reasonably below the frame threshold
Vmax_frac = 0.40;                 % tighten/loosen as needed (0.3–0.5 typical)
Vmax = floor(Vmax_frac * L);

% --- Calibrate tauP (Permutation)
FAR_est = inf;
VmeanP_est = inf;

for tau = tauP_grid
    cfg_tmp = cfg;
    cfg_tmp.tauP = tau;

    acc = 0;
    Vsum = 0;

    for t = 1:MC_cal
        o = run_one_trial_auth(cfg_tmp, nodes, u0, SNR_cal, false); % H0 trials
        acc  = acc  + o.accept_perm;
        Vsum = Vsum + o.Vperm;
    end

    FAR_here = acc/MC_cal;
    Vmean = Vsum/MC_cal;

    if (FAR_here <= targetFAR) && (Vmean <= Vmax)
        best_tauP = tau;
        FAR_est = FAR_here;
        VmeanP_est = Vmean;
        break;  % smallest tau meeting BOTH constraints
    end
end

% --- Calibrate tauIM (IM)
FAR_est_im = inf;
VmeanIM_est = inf;

for tau = tauIM_grid
    cfg_tmp = cfg;
    cfg_tmp.tauIM = tau;

    acc = 0;
    Vsum = 0;

    for t = 1:MC_cal
        o = run_one_trial_auth(cfg_tmp, nodes, u0, SNR_cal, false); % H0 trials
        acc  = acc  + o.accept_im;
        Vsum = Vsum + o.Vim;
    end

    FAR_here = acc/MC_cal;
    Vmean = Vsum/MC_cal;

    if (FAR_here <= targetFAR) && (Vmean <= Vmax)
        best_tauIM = tau;
        FAR_est_im = FAR_here;
        VmeanIM_est = Vmean;
        break;
    end
end

cfg.tauP  = best_tauP;
cfg.tauIM = best_tauIM;

fprintf('Calibrated at SNR=%g dB, target FAR=%.2g, Vmax=%d/%d:\n', SNR_cal, targetFAR, Vmax, L);
fprintf('  tauP=%.4f (FAR~%.3f, mean H0 votes=%.2f)\n', cfg.tauP,  FAR_est,    VmeanP_est);
fprintf('  tauIM=%.4f (FAR~%.3f, mean H0 votes=%.2f)\n', cfg.tauIM, FAR_est_im, VmeanIM_est);

% -------------------------------------------------------------
tauP_list  = cfg.tauP  * [0.85 1.00 1.15];
%tauIM_list = cfg.tauIM * [0.85 1.00 1.15];  % if tauIM is negative, use additive offsets instead
tauIM_list = cfg.tauIM + [-0.05 0 0.05];

FARp = zeros(numel(tauP_list), numel(SNRvec));
FRRp = zeros(numel(tauP_list), numel(SNRvec));
FARi = zeros(numel(tauIM_list), numel(SNRvec));
FRRi = zeros(numel(tauIM_list), numel(SNRvec));

% ---- DEBUG: H0 vote distributions at low vs high SNR ----
dbgSNR = [5 25];
MC_dbg = 500;
Tp_all = [];
Tim_all = [];

for s = dbgSNR
    Vp = zeros(1, MC_dbg);
    Vi = zeros(1, MC_dbg);
    for t = 1:MC_dbg
        o = run_one_trial_auth(cfg, nodes, u0, s, false); % impostor trials only
        Vp(t) = o.Vperm;
        Vi(t) = o.Vim;
        Tp_all  = [Tp_all o.Tp];
        Tim_all = [Tim_all o.Tim];
    end

    need = ceil(cfg.alpha * cfg.Lsym);

    fprintf('\nH0 vote stats at SNRref=%g dB (need >= %d of %d):\n', s, need, cfg.Lsym);
    fprintf(' Perm votes: mean=%.2f, min=%d, max=%d, P(V>=need)=%.3f\n', ...
        mean(Vp), min(Vp), max(Vp), mean(Vp >= need));
    fprintf(' IM   votes: mean=%.2f, min=%d, max=%d, P(V>=need)=%.3f\n', ...
        mean(Vi), min(Vi), max(Vi), mean(Vi >= need));

    fprintf(' H0 per-symbol |Tperm|: mean=%.4f std=%.4f\n', mean(Tp_all), std(Tp_all));
    fprintf(' H0 per-symbol Tim:     mean=%.4f std=%.4f\n', mean(Tim_all), std(Tim_all));

    % Optional: show top few counts to see if it clusters at a value
    tabP = tabulate(Vp(:));
    tabI = tabulate(Vi(:));
    tabP = sortrows(tabP, -3); tabI = sortrows(tabI, -3);
    fprintf(' Perm most common V counts (value, count, pct):\n');
    disp(tabP(1:min(5,size(tabP,1)),:));
    fprintf(' IM most common V counts (value, count, pct):\n');
    disp(tabI(1:min(5,size(tabI,1)),:));
end
% ---------------------------------------------------------


%FARp = zeros(size(SNRvec)); FRRp = zeros(size(SNRvec));
%FARi = zeros(size(SNRvec)); FRRi = zeros(size(SNRvec));

% Allocate result arrays (multi-tau)

for si = 1:numel(SNRvec)
    snr = SNRvec(si);

    % ---- Permutation curves for each tauP
    for ti = 1:numel(tauP_list)
        cfg_tmp = cfg;
        cfg_tmp.tauP = tauP_list(ti);

        % legit -> FRR
        rej = 0;
        for t = 1:MC
            o = run_one_trial_auth(cfg_tmp, nodes, u0, snr, true);
            rej = rej + (~o.accept_perm);
        end
        FRRp(ti, si) = rej/MC;

        % impostor -> FAR
        acc = 0;
        for t = 1:MC
            o = run_one_trial_auth(cfg_tmp, nodes, u0, snr, false);
            acc = acc + (o.accept_perm);
        end
        FARp(ti, si) = acc/MC;
    end

    % ---- IM curves for each tauIM
    for ti = 1:numel(tauIM_list)
        cfg_tmp = cfg;
        cfg_tmp.tauIM = tauIM_list(ti);

        % legit -> FRR
        rej = 0;
        for t = 1:MC
            o = run_one_trial_auth(cfg_tmp, nodes, u0, snr, true);
            rej = rej + (~o.accept_im);
        end
        FRRi(ti, si) = rej/MC;

        % impostor -> FAR
        acc = 0;
        for t = 1:MC
            o = run_one_trial_auth(cfg_tmp, nodes, u0, snr, false);
            acc = acc + (o.accept_im);
        end
        FARi(ti, si) = acc/MC;
    end

    % Print just the calibrated tau (middle entry) for quick logging
    midP  = ceil(numel(tauP_list)/2);
    midIM = ceil(numel(tauIM_list)/2);
    fprintf('SNRref=%g dB: FRR(P)=%.3f FAR(P)=%.3f | FRR(IM)=%.3f FAR(IM)=%.3f\n', ...
        snr, FRRp(midP,si), FARp(midP,si), FRRi(midIM,si), FARi(midIM,si));
end


% Plot (simple)
% figure; plot(SNRvec, FARp, '-o'); hold on; plot(SNRvec, FRRp, '-s');
% xlabel('SNR_{ref} at d_0 (dB)'); ylabel('Probability'); grid on;
% legend('FAR (Perm)','FRR (Perm)','Location','best');
% title('Permutation Authentication: FAR/FRR vs SNR');
% 
% figure; plot(SNRvec, FARi, '-o'); hold on; plot(SNRvec, FRRi, '-s');
% xlabel('SNR_{ref} at d_0 (dB)'); ylabel('Probability'); grid on;
% legend('FAR (IM)','FRR (IM)','Location','best');
% title('IM Authentication: FAR/FRR vs SNR');

figure; hold on; grid on;
for ti = 1:numel(tauP_list)
    plot(SNRvec, FARp(ti,:), '-o');
    plot(SNRvec, FRRp(ti,:), '--s');
end
xlabel('SNR_{ref} at d_0 (dB)'); ylabel('Probability');
title('Permutation Authentication: FAR/FRR vs SNR for multiple \tau_P');

% Build legend entries
leg = strings(1, 2*numel(tauP_list));
for ti = 1:numel(tauP_list)
    leg(2*ti-1) = sprintf('FAR (\\tau_P=%.3f)', tauP_list(ti));
    leg(2*ti)   = sprintf('FRR (\\tau_P=%.3f)', tauP_list(ti));
end
legend(leg, 'Location','best');

figure; hold on; grid on;
for ti = 1:numel(tauIM_list)
    plot(SNRvec, FARi(ti,:), '-o');
    plot(SNRvec, FRRi(ti,:), '--s');
end
xlabel('SNR_{ref} at d_0 (dB)'); ylabel('Probability');
title('IM Authentication: FAR/FRR vs SNR for multiple \tau_{IM}');

leg = strings(1, 2*numel(tauIM_list));
for ti = 1:numel(tauIM_list)
    leg(2*ti-1) = sprintf('FAR (\\tau_{IM}=%.3f)', tauIM_list(ti));
    leg(2*ti)   = sprintf('FRR (\\tau_{IM}=%.3f)', tauIM_list(ti));
end
legend(leg, 'Location','best');

end
