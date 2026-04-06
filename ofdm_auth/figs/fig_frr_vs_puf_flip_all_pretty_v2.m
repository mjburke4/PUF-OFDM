function fig_frr_vs_puf_flip_all_pretty_v2()

%% =========================
% Common settings
% ==========================
alpha = 0.75;
L_list = [4 64];
SNR_list = [0 15];
MC = 5000;
MC_cal = 3000;

p_list = [0 0.002 0.005 0.01 0.02 0.05];

% Calibration settings
SNR_cal = 15;          % reference SNR for fixed-threshold calibration
L_cal   = 8;           % reference L for fixed-threshold calibration
p_cal   = 0;           % calibrate using stable PUF first
targetFAR = 0.01;      % security target

U = 10;
m = 64;
u0 = 3;

nodes = init_node_population(U, m, 'seed', 7);

% Fix distance so SNRref is meaningful for perm/IM branch
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

% --- Calibration for Perm / IM
cfg_cal = cfg;
cfg_cal.Lsym = L_cal;

tauP_grid  = linspace(0.00, 0.50, 101);
tauIM_grid = linspace(0.00, 1.00, 121);

% Calibrate tauP: among thresholds with FAR <= target, choose min FRR
best_tauP = tauP_grid(end);
best_FRR_P = inf;
best_FAR_P = inf;

for tau = tauP_grid
    cfgTmp = cfg_cal;
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
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, true, p_cal);
        rej = rej + (~o.accept_perm);
    end
    FRR_here = rej / MC_cal;

    if (FAR_here <= targetFAR) && (FRR_here < best_FRR_P)
        best_tauP = tau;
        best_FRR_P = FRR_here;
        best_FAR_P = FAR_here;
    end
end

% Fallback if no threshold satisfies target FAR
if isinf(best_FRR_P)
    warning('No tauP met FAR target; using approximate EER-like choice.');
    best_gap = inf;
    for tau = tauP_grid
        cfgTmp = cfg_cal;
        cfgTmp.tauP = tau;

        acc = 0; rej = 0;
        for t = 1:MC_cal
            o0 = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, false, 0);
            o1 = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, true, p_cal);
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

% Calibrate tauIM: among thresholds with FAR <= target, choose min FRR
best_tauIM = tauIM_grid(end);
best_FRR_IM = inf;
best_FAR_IM = inf;

for tau = tauIM_grid
    cfgTmp = cfg_cal;
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
        o = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, true, p_cal);
        rej = rej + (~o.accept_im);
    end
    FRR_here = rej / MC_cal;

    if (FAR_here <= targetFAR) && (FRR_here < best_FRR_IM)
        best_tauIM = tau;
        best_FRR_IM = FRR_here;
        best_FAR_IM = FAR_here;
    end
end

% Fallback if no threshold satisfies target FAR
if isinf(best_FRR_IM)
    warning('No tauIM met FAR target; using approximate EER-like choice.');
    best_gap = inf;
    for tau = tauIM_grid
        cfgTmp = cfg_cal;
        cfgTmp.tauIM = tau;

        acc = 0; rej = 0;
        for t = 1:MC_cal
            o0 = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, false, 0);
            o1 = run_one_trial_auth(cfgTmp, nodes, u0, SNR_cal, true, p_cal);
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

fprintf('\nCalibration at SNR=%g dB, L=%d, p_flip=%g:\n', SNR_cal, L_cal, p_cal);
fprintf('  tauP  = %.4f   (FAR=%.4f, FRR=%.4f)\n', cfg.tauP,  best_FAR_P,  best_FRR_P);
fprintf('  tauIM = %.4f   (FAR=%.4f, FRR=%.4f)\n', cfg.tauIM, best_FAR_IM, best_FRR_IM);

%% =========================
% Sweep Perm / IM
% ==========================
FRRp = zeros(numel(L_list), numel(SNR_list), numel(p_list));
FRRi = zeros(numel(L_list), numel(SNR_list), numel(p_list));

for li = 1:numel(L_list)
    cfgL = cfg;
    cfgL.Lsym = L_list(li);

    for si = 1:numel(SNR_list)
        snr = SNR_list(si);

        for pi = 1:numel(p_list)
            p_flip = p_list(pi);

            rejP = 0;
            rejI = 0;

            for t = 1:MC
                o = run_one_trial_auth(cfgL, nodes, u0, snr, true, p_flip);
                rejP = rejP + (~o.accept_perm);
                rejI = rejI + (~o.accept_im);
            end

            FRRp(li, si, pi) = rejP / MC;
            FRRi(li, si, pi) = rejI / MC;
        end
    end
end

%% =========================
% DPPA config
% ==========================
cfgD = cfg_default_dppa();
cfgD.chan.type = 'awgn';
cfgD.alpha = alpha;

% Optional for WLAN-like consistency:
% cfgD.Npilots = 4;

cfgD_cal = cfgD;
cfgD_cal.Lsym = L_cal;

tauD_grid = linspace(0.60, 1.00, 81);

best_tauD = tauD_grid(end);
best_FRR_D = inf;
best_FAR_D = inf;

for tau = tauD_grid
    cfgTmp = cfgD_cal;
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
        o = run_one_trial_dppa(cfgTmp, nodes, u0, SNR_cal, true, p_cal);
        rej = rej + (~o.accept_dppa);
    end
    FRR_here = rej / MC_cal;

    if (FAR_here <= targetFAR) && (FRR_here < best_FRR_D)
        best_tauD = tau;
        best_FRR_D = FRR_here;
        best_FAR_D = FAR_here;
    end
end

% Fallback if no threshold satisfies target FAR
if isinf(best_FRR_D)
    warning('No tauD met FAR target; using approximate EER-like choice.');
    best_gap = inf;
    for tau = tauD_grid
        cfgTmp = cfgD_cal;
        cfgTmp.tauD = tau;

        acc = 0; rej = 0;
        for t = 1:MC_cal
            o0 = run_one_trial_dppa(cfgTmp, nodes, u0, SNR_cal, false, 0);
            o1 = run_one_trial_dppa(cfgTmp, nodes, u0, SNR_cal, true, p_cal);
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

fprintf('  tauD  = %.4f   (FAR=%.4f, FRR=%.4f)\n\n', cfgD.tauD, best_FAR_D, best_FRR_D);

%% =========================
% Sweep DPPA
% ==========================
FRRd = zeros(numel(L_list), numel(SNR_list), numel(p_list));

for li = 1:numel(L_list)
    cfgDL = cfgD;
    cfgDL.Lsym = L_list(li);

    for si = 1:numel(SNR_list)
        snr = SNR_list(si);

        for pi = 1:numel(p_list)
            p_flip = p_list(pi);

            rejD = 0;
            for t = 1:MC
                o = run_one_trial_dppa(cfgDL, nodes, u0, snr, true, p_flip);
                rejD = rejD + (~o.accept_dppa);
            end

            FRRd(li, si, pi) = rejD / MC;
        end
    end
end

%% =========================
% Pretty plot
% ==========================
lw = 1.9;
ms = 7;
axFS = 12;
titleFS = 13;
sgFS = 14;

fig = figure('Color','w', 'Position', [100 100 1350 420]);
tl = tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');

methodNames = {'Permutation', 'IM', 'DPPA'};
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

    for li = 1:numel(L_list)
        for si = 1:numel(SNR_list)
            h = plot(ax, p_list, max(squeeze(allFRR{k}(li,si,:))', 1/MC), ...
                '-o', 'LineWidth', lw, 'MarkerSize', ms);

            if k == 1
                legendHandles = [legendHandles h]; %#ok<AGROW>
            end
        end
    end

    xlabel(ax, 'PUF bit-flip probability, p_{flip}', 'FontSize', axFS);
    if k == 1
        ylabel(ax, 'FRR', 'FontSize', axFS);
    end

    title(ax, sprintf('(%c) %s', 'a'+(k-1), methodNames{k}), ...
        'FontWeight', 'bold', 'FontSize', titleFS);

    xlim(ax, [min(p_list) max(p_list)]);
    xticks(ax, p_list);
    ylim(ax, [1/MC 1]);
end

sgtitle(tl, sprintf('FRR vs PUF Bit-Flip Probability (\\alpha = %.2f, calibrated at SNR = %d dB, L = %d)', ...
    alpha, SNR_cal, L_cal), ...
    'FontWeight', 'bold', 'FontSize', sgFS);

legendStrings = {};
for li = 1:numel(L_list)
    for si = 1:numel(SNR_list)
        legendStrings{end+1} = sprintf('L = %d, SNR = %d dB', L_list(li), SNR_list(si));
    end
end

lgd = legend(legendHandles, legendStrings, ...
    'Orientation', 'horizontal', ...
    'Location', 'southoutside', ...
    'NumColumns', numel(legendStrings), ...
    'FontSize', 11);
lgd.Layout.Tile = 'south';

end