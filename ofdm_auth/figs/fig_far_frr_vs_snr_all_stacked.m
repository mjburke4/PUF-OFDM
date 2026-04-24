function fig_far_frr_vs_snr_all_stacked()
addpath 'C:\Users\michael.j.burke10\OneDrive - US Army\Desktop\Grad school stuff\PUF based dynamic spreading code generation\simulations\PUF_Modulation_ofdm_Sim\PUF-OFDM\ofdm_auth\lib'
addpath 'C:\Users\michael.j.burke10\OneDrive - US Army\Desktop\Grad school stuff\PUF based dynamic spreading code generation\simulations\PUF_Modulation_ofdm_Sim\PUF-OFDM\ofdm_auth\cfg'
addpath 'C:\Users\michael.j.burke10\OneDrive - US Army\Desktop\Grad school stuff\PUF based dynamic spreading code generation\simulations\PUF_Modulation_ofdm_Sim\PUF-OFDM\ofdm_auth\sim'
%% =========================
% Common settings
% ==========================
alpha = 0.75;
L_cal = 16;              % reference frame length for calibration
SNR_cal = 15;           % calibration SNR
targetFAR = 0.01;       % fixed-threshold operating target
MC = 5000;              % per-SNR trials
MC_cal = 3000;          % calibration trials

SNRvec = [-5 5 15 25];

U = 10;
m = 64;
u0 = 3;

nodes = init_node_population(U, m, 'seed', 7);

% Fix distance so SNRref is meaningful
d0 = 10;
for u = 1:U
    nodes(u).distance_m = d0;
end

%% =========================
% Perm / IM config
% ==========================
cfg = cfg_default();
cfg.chan.type = 'awgn';
cfg.alpha = alpha;
cfg.Lsym = L_cal;

tauP_grid  = linspace(0.00, 0.50, 101);
tauIM_grid = linspace(0.00, 1.00, 121);

% -------- Calibrate tauP --------
best_tauP = tauP_grid(end);
best_FRR_P = inf;
best_FAR_P = inf;

for tau = tauP_grid
    cfgTmp = cfg;
    cfgTmp.tauP = tau;

    % FAR on H0
    acc = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, false, 0);
        acc = acc + o.accept_perm;
    end
    FAR_here = acc / MC_cal;

    % FRR on H1
    rej = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, true, 0);
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
    warning('No tauP met FAR target; using approximate EER fallback.');
    best_gap = inf;
    for tau = tauP_grid
        cfgTmp = cfg;
        cfgTmp.tauP = tau;

        acc = 0; rej = 0;
        for t = 1:MC_cal
            o0 = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, false, 0);
            o1 = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, true, 0);
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

% -------- Calibrate tauIM --------
best_tauIM = tauIM_grid(end);
best_FRR_IM = inf;
best_FAR_IM = inf;

for tau = tauIM_grid
    cfgTmp = cfg;
    cfgTmp.tauIM = tau;

    % FAR on H0
    acc = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, false, 0);
        acc = acc + o.accept_im;
    end
    FAR_here = acc / MC_cal;

    % FRR on H1
    rej = 0;
    for t = 1:MC_cal
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, true, 0);
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
    warning('No tauIM met FAR target; using approximate EER fallback.');
    best_gap = inf;
    for tau = tauIM_grid
        cfgTmp = cfg;
        cfgTmp.tauIM = tau;

        acc = 0; rej = 0;
        for t = 1:MC_cal
            o0 = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, false, 0);
            o1 = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, true, 0);
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

cfg.tauP = best_tauP;
cfg.tauIM = best_tauIM;

%% =========================
% DPPA config
% ==========================
cfgD = cfg_default_dppa();
cfgD.chan.type = 'awgn';
cfgD.alpha = alpha;
cfgD.Lsym = L_cal;

% Easy switch if you want to test richer alphabet:
% cfgD.phi_choices_deg = [-20 -10 10 20];

tauD_grid = linspace(0.50, 0.999, 120);

best_tauD = tauD_grid(end);
best_FRR_D = inf;
best_FAR_D = inf;

for tau = tauD_grid
    cfgTmp = cfgD;
    cfgTmp.tauD = tau;

    % FAR on H0
    acc = 0;
    for t = 1:MC_cal
        o = run_one_trial_dppa(cfgTmp, nodes, u0, SNR_cal, false, 0);
        acc = acc + o.accept_dppa;
    end
    FAR_here = acc / MC_cal;

    % FRR on H1
    rej = 0;
    for t = 1:MC_cal
        o = run_one_trial_dppa(cfgTmp, nodes, u0, SNR_cal, true, 0);
        rej = rej + (~o.accept_dppa);
    end
    FRR_here = rej / MC_cal;

    if (FAR_here <= targetFAR) && (FRR_here < best_FRR_D)
        best_tauD = tau;
        best_FRR_D = FRR_here;
        best_FAR_D = FAR_here;
    end
end

if isinf(best_FRR_D)
    warning('No tauD met FAR target; using approximate EER fallback.');
    best_gap = inf;
    for tau = tauD_grid
        cfgTmp = cfgD;
        cfgTmp.tauD = tau;

        acc = 0; rej = 0;
        for t = 1:MC_cal
            o0 = run_one_trial_dppa(cfgTmp, nodes, u0, SNR_cal, false, 0);
            o1 = run_one_trial_dppa(cfgTmp, nodes, u0, SNR_cal, true, 0);
            acc = acc + o0.accept_dppa;
            rej = rej + (~o1.accept_dppa);
        end
        FAR_here = acc / MC_cal;
        FRR_here = rej / MC_cal;
        gap = abs(FAR_here - FRR_here);

        if gap < best_gap
            best_gap = gap;
            best_tauD = tau;
            best_FRR_D = FRR_here;
            best_FAR_D = FAR_here;
        end
    end
end

cfgD.tauD = best_tauD;

fprintf('\nCalibrated thresholds at SNR=%g dB, L=%d:\n', SNR_cal, L_cal);
fprintf('  tauP  = %.4f   (FAR=%.4f, FRR=%.4f)\n', cfg.tauP, best_FAR_P, best_FRR_P);
fprintf('  tauIM = %.4f   (FAR=%.4f, FRR=%.4f)\n', cfg.tauIM, best_FAR_IM, best_FRR_IM);
fprintf('  tauD  = %.4f   (FAR=%.4f, FRR=%.4f)\n\n', cfgD.tauD, best_FAR_D, best_FRR_D);

%% =========================
% Evaluate across SNR
% ==========================
FARp = zeros(size(SNRvec));  FRRp = zeros(size(SNRvec));
FARi = zeros(size(SNRvec));  FRRi = zeros(size(SNRvec));
FARd = zeros(size(SNRvec));  FRRd = zeros(size(SNRvec));

for si = 1:numel(SNRvec)
    snr = SNRvec(si);

    % -------- Permutation
    rej = 0; acc = 0;
    for t = 1:MC
        o1 = run_one_trial_auth(cfg, nodes, u0, snr, true, 0);
        o0 = run_one_trial_auth(cfg, nodes, u0, snr, false, 0);
        rej = rej + (~o1.accept_perm);
        acc = acc + o0.accept_perm;
    end
    FRRp(si) = rej / MC;
    FARp(si) = acc / MC;

    % -------- IM
    rej = 0; acc = 0;
    for t = 1:MC
        o1 = run_one_trial_auth(cfg, nodes, u0, snr, true, 0);
        o0 = run_one_trial_auth(cfg, nodes, u0, snr, false, 0);
        rej = rej + (~o1.accept_im);
        acc = acc + o0.accept_im;
    end
    FRRi(si) = rej / MC;
    FARi(si) = acc / MC;

    % -------- DPPA
    rej = 0; acc = 0;
    for t = 1:MC
        o1 = run_one_trial_dppa(cfgD, nodes, u0, snr, true, 0);
        o0 = run_one_trial_dppa(cfgD, nodes, u0, snr, false, 0);
        rej = rej + (~o1.accept_dppa);
        acc = acc + o0.accept_dppa;
    end
    FRRd(si) = rej / MC;
    FARd(si) = acc / MC;

    fprintf(['SNR=%g dB: ' ...
             'Perm(FAR=%.4f, FRR=%.4f) | ' ...
             'IM(FAR=%.4f, FRR=%.4f) | ' ...
             'DPPA(FAR=%.4f, FRR=%.4f)\n'], ...
             snr, FARp(si), FRRp(si), FARi(si), FRRi(si), FARd(si), FRRd(si));
end

%% =========================
% Pretty stacked plot
% ==========================
lw = 2.0;
ms = 7;
axFS = 12;
titleFS = 13;
sgFS = 14;

fig = figure('Color','w', 'Position', [100 60 900 980]);
tl = tiledlayout(3,1, 'TileSpacing','compact', 'Padding','compact');

methodNames = {'Permutation', 'Index Modulation', 'DPPA'};
allFAR = {FARp, FARi, FARd};
allFRR = {FRRp, FRRi, FRRd};

legendHandles = [];

for k = 1:3
    ax = nexttile;
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');

    set(ax, 'YScale', 'log', ...
        'FontSize', axFS, ...
        'LineWidth', 1.0, ...
        'XMinorGrid', 'off', ...
        'YMinorGrid', 'on');

    h1 = plot(ax, SNRvec, max(allFAR{k}, 1/MC), '-o', ...
        'LineWidth', lw, 'MarkerSize', ms);
    h2 = plot(ax, SNRvec, max(allFRR{k}, 1/MC), '--s', ...
        'LineWidth', lw, 'MarkerSize', ms);

    if k == 1
        legendHandles = [h1 h2];
    end

    xlabel(ax, 'SNR_{ref} (dB)', 'FontSize', axFS);
    ylabel(ax, 'Probability', 'FontSize', axFS);
    title(ax, sprintf('(%c) %s', 'a'+(k-1), methodNames{k}), ...
        'FontWeight', 'bold', 'FontSize', titleFS);

    xticks(ax, SNRvec);
    xlim(ax, [min(SNRvec) max(SNRvec)]);
    ylim(ax, [1/MC 1]);
end

sgtitle(tl, sprintf('FAR/FRR vs SNR (\\alpha = %.2f, calibrated at SNR = %d dB, L = %d)', ...
    alpha, SNR_cal, L_cal), ...
    'FontWeight', 'bold', 'FontSize', sgFS);

lgd = legend(legendHandles, {'FAR', 'FRR'}, ...
    'Orientation', 'horizontal', ...
    'Location', 'southoutside', ...
    'NumColumns', 2, ...
    'FontSize', 11);
lgd.Layout.Tile = 'south';

end