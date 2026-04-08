function fig_far_frr_vs_L_all_pretty_v2()

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

% Plot styling
lw = 1.8;
ms = 7;
axFS = 12;
titleFS = 13;
sgFS = 14;

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

Vmax = floor(Vmax_frac * cfg_cal.Lsym);

tauP_grid  = linspace(0.00, 0.50, 101);
tauIM_grid = linspace(0.00, 0.50, 101);

% --- Calibrate tauP
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

% --- Calibrate tauIM
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

% For paper consistency, consider setting this to 4 if desired:
cfgD.Npilots = 4;

FARd = zeros(numel(SNR_list), numel(L_list));
FRRd = zeros(numel(SNR_list), numel(L_list));

cfgD_cal = cfgD;
cfgD_cal.Lsym = L_cal;

Ldiff_cal = cfgD_cal.Lsym - 1;
VmaxD     = floor(Vmax_frac * Ldiff_cal);

tauD_grid = linspace(0.60, 1.00, 81);
best_tauD   = tauD_grid(1);
bestFRR_cal = inf;
bestFAR_cal = inf;

for tau = tauD_grid
    cfgTmp = cfgD_cal;
    cfgTmp.tauD = tau;

    % H0 -> FAR
    acc = 0;
    for t = 1:MC_cal
        o = run_one_trial_dppa(cfgTmp, nodes, u0, SNR_cal, false);
        acc = acc + o.accept_dppa;
    end
    FAR_here = acc / MC_cal;

    % H1 -> FRR
    rej = 0;
    for t = 1:MC_cal
        o = run_one_trial_dppa(cfgTmp, nodes, u0, SNR_cal, true);
        rej = rej + (~o.accept_dppa);
    end
    FRR_here = rej / MC_cal;

    % choose the lowest FRR among thresholds that satisfy FAR target
    if (FAR_here <= targetFAR) && (FRR_here < bestFRR_cal)
        best_tauD   = tau;
        bestFRR_cal = FRR_here;
        bestFAR_cal = FAR_here;
    end
end

fprintf('DPPA calibrated tauD=%.4f, FAR=%.4f, FRR=%.4f at %g dB\n', ...
    best_tauD, bestFAR_cal, bestFRR_cal, SNR_cal);
% best_tauD = tauD_grid(end);
% for tau = tauD_grid
%     cfgTmp = cfgD_cal;
%     cfgTmp.tauD = tau;
% 
%     acc = 0; Vsum = 0;
%     for t = 1:MC_cal
%         o = run_one_trial_dppa(cfgTmp, nodes, u0, SNR_cal, false);
%         acc  = acc  + o.accept_dppa;
%         Vsum = Vsum + o.Vd_H1;
%     end
%     FAR_here = acc/MC_cal;
%     Vmean    = Vsum/MC_cal;
% 
%     if (FAR_here <= targetFAR) && (Vmean <= VmaxD)
%         best_tauD = tau;
%         break;
%     end
% end

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
% Plot
% ==========================
fig = figure('Color','w', 'Position', [100 100 1350 420]);
tl = tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');

methodNames = {'Permutation','IM','DPPA'};
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

    % Plot curves
    for si = 1:numel(SNR_list)
        h1 = plot(ax, L_list, max(allFAR{k}(si,:), 1/MC), '-o', ...
            'LineWidth', lw, 'MarkerSize', ms);
        h2 = plot(ax, L_list, max(allFRR{k}(si,:), 1/MC), '--s', ...
            'LineWidth', lw, 'MarkerSize', ms);

        if k == 1
            legendHandles = [legendHandles h1 h2]; %#ok<AGROW>
        end
    end

    xlabel(ax, 'Frame length, L', 'FontSize', axFS);
    if k == 1
        ylabel(ax, 'Probability', 'FontSize', axFS);
    end

    title(ax, sprintf('(%c) %s', 'a'+(k-1), methodNames{k}), ...
        'FontWeight', 'bold', 'FontSize', titleFS);

    xlim(ax, [min(L_list) max(L_list)]);
    xticks(ax, L_list);

    % Optional y-limits if you want more uniform appearance
    ylim(ax, [1/MC 1]);
end

sgtitle(tl, sprintf('FAR/FRR vs Frame Length (\\alpha = %.2f, calibrated at SNR = %d dB, L = %d)', ...
    alpha, SNR_cal, L_cal), ...
    'FontWeight', 'bold', 'FontSize', sgFS);

legendStrings = {};
for si = 1:numel(SNR_list)
    legendStrings{end+1} = sprintf('FAR (%d dB)', SNR_list(si));
    legendStrings{end+1} = sprintf('FRR (%d dB)', SNR_list(si));
end

lgd = legend(legendHandles, legendStrings, ...
    'Orientation', 'horizontal', ...
    'Location', 'southoutside', ...
    'NumColumns', 3, ...
    'FontSize', 11);
lgd.Layout.Tile = 'south';

end