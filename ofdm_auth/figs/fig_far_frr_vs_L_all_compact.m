function fig_far_frr_vs_L_all_compact()

%% =========================
% Common sweep settings
% ==========================
SNR_list = [-5 5 15];
L_list   = [4 8 16 32 64];
MC       = 6000;

alpha     = 0.75;
L_cal     = 8;
SNR_cal   = 15;
MC_cal    = 3000;
targetFAR = 0.01;
Vmax_frac = 0.40;

U = 10; 
m = 64;
nodes = init_node_population(U, m, 'seed', 7);
u0 = 3;

%% =========================
% PERM / IM
% ==========================
cfg = cfg_default();
cfg.chan.type = 'awgn';
cfg.alpha = alpha;

d0 = 10;
for u = 1:U
    nodes(u).distance_m = d0;
end

FARp = zeros(numel(SNR_list), numel(L_list));
FRRp = zeros(numel(SNR_list), numel(L_list));
FARi = zeros(numel(SNR_list), numel(L_list));
FRRi = zeros(numel(SNR_list), numel(L_list));

cfg_cal = cfg;
cfg_cal.Lsym = L_cal;

need_cal = ceil(cfg_cal.alpha * cfg_cal.Lsym);
Vmax     = floor(Vmax_frac * cfg_cal.Lsym);

tauP_grid  = linspace(0.00, 0.50, 101);
tauIM_grid = linspace(0.00, 0.50, 101);

% Calibrate tauP
best_tauP = tauP_grid(end);
for tau = tauP_grid
    cfgTmp = cfg_cal;
    cfgTmp.tauP = tau;

    acc = 0; Vsum = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, false);
        acc  = acc  + o.accept_perm;
        Vsum = Vsum + o.Vperm;
    end
    FAR_here = acc/MC_cal;
    Vmean    = Vsum/MC_cal;

    if (FAR_here <= targetFAR) && (Vmean <= Vmax)
        best_tauP = tau;
        break;
    end
end

% Calibrate tauIM
best_tauIM = tauIM_grid(end);
for tau = tauIM_grid
    cfgTmp = cfg_cal;
    cfgTmp.tauIM = tau;

    acc = 0; Vsum = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, false);
        acc  = acc  + o.accept_im;
        Vsum = Vsum + o.Vim;
    end
    FAR_here = acc/MC_cal;
    Vmean    = Vsum/MC_cal;

    if (FAR_here <= targetFAR) && (Vmean <= Vmax)
        best_tauIM = tau;
        break;
    end
end

cfg.tauP  = best_tauP;
cfg.tauIM = best_tauIM;

for li = 1:numel(L_list)
    cfgL = cfg;
    cfgL.Lsym = L_list(li);

    for si = 1:numel(SNR_list)
        snr = SNR_list(si);

        % Perm FAR
        acc = 0;
        for t = 1:MC
            o = run_one_trial_auth(cfgL, nodes, u0, snr, false);
            acc = acc + o.accept_perm;
        end
        FARp(si,li) = acc/MC;

        % Perm FRR
        rej = 0;
        for t = 1:MC
            o = run_one_trial_auth(cfgL, nodes, u0, snr, true);
            rej = rej + (~o.accept_perm);
        end
        FRRp(si,li) = rej/MC;

        % IM FAR
        acc = 0;
        for t = 1:MC
            o = run_one_trial_auth(cfgL, nodes, u0, snr, false);
            acc = acc + o.accept_im;
        end
        FARi(si,li) = acc/MC;

        % IM FRR
        rej = 0;
        for t = 1:MC
            o = run_one_trial_auth(cfgL, nodes, u0, snr, true);
            rej = rej + (~o.accept_im);
        end
        FRRi(si,li) = rej/MC;
    end
end

%% =========================
% DPPA
% ==========================
cfgD = cfg_default_dppa();
cfgD.chan.type = 'awgn';
cfgD.alpha = alpha;

FARd = zeros(numel(SNR_list), numel(L_list));
FRRd = zeros(numel(SNR_list), numel(L_list));

cfgD_cal = cfgD;
cfgD_cal.Lsym = L_cal;

Ldiff_cal = cfgD_cal.Lsym - 1;
VmaxD     = floor(Vmax_frac * Ldiff_cal);

tauD_grid = linspace(0.60, 1.00, 81);

best_tauD = tauD_grid(end);
for tau = tauD_grid
    cfgTmp = cfgD_cal;
    cfgTmp.tauD = tau;

    acc = 0; Vsum = 0;
    for t = 1:MC_cal
        o = run_one_trial_dppa(cfgTmp, nodes, u0, SNR_cal, false);
        acc  = acc  + o.accept_dppa;
        Vsum = Vsum + o.Vd_H1;
    end
    FAR_here = acc/MC_cal;
    Vmean    = Vsum/MC_cal;

    if (FAR_here <= targetFAR) && (Vmean <= VmaxD)
        best_tauD = tau;
        break;
    end
end

cfgD.tauD = best_tauD;

for li = 1:numel(L_list)
    cfgL = cfgD;
    cfgL.Lsym = L_list(li);

    for si = 1:numel(SNR_list)
        snr = SNR_list(si);

        % DPPA FAR
        acc = 0;
        for t = 1:MC
            o = run_one_trial_dppa(cfgL, nodes, u0, snr, false);
            acc = acc + o.accept_dppa;
        end
        FARd(si,li) = acc/MC;

        % DPPA FRR
        rej = 0;
        for t = 1:MC
            o = run_one_trial_dppa(cfgL, nodes, u0, snr, true);
            rej = rej + (~o.accept_dppa);
        end
        FRRd(si,li) = rej/MC;
    end
end

%% =========================
% Compact combined plot
% ==========================
figure;
tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');

% -------- Permutation
nexttile; hold on; grid on; set(gca,'YScale','log');
for si = 1:numel(SNR_list)
    plot(L_list, max(FARp(si,:), 1/MC), '-o', 'LineWidth', 1.8);
    plot(L_list, max(FRRp(si,:), 1/MC), '--s', 'LineWidth', 1.8);
end
xlabel('L');
ylabel('Probability');
title('Permutation');

% -------- IM
nexttile; hold on; grid on; set(gca,'YScale','log');
for si = 1:numel(SNR_list)
    plot(L_list, max(FARi(si,:), 1/MC), '-o', 'LineWidth', 1.8);
    plot(L_list, max(FRRi(si,:), 1/MC), '--s', 'LineWidth', 1.8);
end
xlabel('L');
ylabel('Probability');
title('IM');

% -------- DPPA
nexttile; hold on; grid on; set(gca,'YScale','log');
for si = 1:numel(SNR_list)
    plot(L_list, max(FARd(si,:), 1/MC), '-o', 'LineWidth', 1.8);
    plot(L_list, max(FRRd(si,:), 1/MC), '--s', 'LineWidth', 1.8);
end
xlabel('L');
ylabel('Probability');
title('DPPA');

sgtitle(sprintf('FAR/FRR vs L (\\alpha=%.2f, MC=%d, fixed thresholds calibrated @ SNR=%gdB, L=%d)', ...
    alpha, MC, SNR_cal, L_cal));

legendStrings = {};
for si = 1:numel(SNR_list)
    legendStrings{end+1} = sprintf('FAR (%g dB)', SNR_list(si));
    legendStrings{end+1} = sprintf('FRR (%g dB)', SNR_list(si));
end
legend(legendStrings, 'Location', 'bestoutside');

end