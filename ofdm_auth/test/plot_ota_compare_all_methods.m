function plot_ota_compare_all_methods()

% Loads OTA result structs from:
%   results_perm.mat
%   results_im.mat
%   results_dppa.mat
%
% Each file should contain a variable named "results".

clear; clc; close all;

files = {
    'Permutation', 'results_perm.mat', 'TpH1_mean',  'TpH0_mean';
    'IM',          'results_im.mat',   'TimH1_mean', 'TimH0_mean';
    'DPPA',        'results_dppa.mat', 'TdH1_mean',  'TdH0_mean';
};

nM = size(files,1);

method = strings(nM,1);
BER_mean = zeros(nM,1);
BER_std  = zeros(nM,1);
H1_mean  = zeros(nM,1);
H0_mean  = zeros(nM,1);
H1_std   = zeros(nM,1);
H0_std   = zeros(nM,1);
gap_mean = zeros(nM,1);
gap_std  = zeros(nM,1);
accH1    = zeros(nM,1);
accH0    = zeros(nM,1);

allBER = cell(nM,1);
allH1  = cell(nM,1);
allH0  = cell(nM,1);
allGap = cell(nM,1);

for k = 1:nM
    method(k) = files{k,1};
    fname     = files{k,2};
    h1field   = files{k,3};
    h0field   = files{k,4};

    S = load(fname);

    if isfield(S,'results')
        R = S.results;
    else
        error('%s does not contain variable named results.', fname);
    end

    % Extract BER
    ber = get_field_vector(R, 'ber');

    % Extract H1/H0 mean scores
    h1 = get_field_vector(R, h1field);
    h0 = get_field_vector(R, h0field);

    % If mean fields are missing, try raw arrays and average per trial
    if isempty(h1)
        h1 = get_raw_metric_fallback(R, method(k), true);
    end
    if isempty(h0)
        h0 = get_raw_metric_fallback(R, method(k), false);
    end

    % Extract accept flags if available
    a1 = get_field_vector(R, 'acceptH1');
    a0 = get_field_vector(R, 'acceptH0');

    gap = h1 - h0;

    allBER{k} = ber;
    allH1{k}  = h1;
    allH0{k}  = h0;
    allGap{k} = gap;

    BER_mean(k) = mean(ber, 'omitnan');
    BER_std(k)  = std(ber,  'omitnan');

    H1_mean(k) = mean(h1, 'omitnan');
    H0_mean(k) = mean(h0, 'omitnan');
    H1_std(k)  = std(h1,  'omitnan');
    H0_std(k)  = std(h0,  'omitnan');

    gap_mean(k) = mean(gap, 'omitnan');
    gap_std(k)  = std(gap,  'omitnan');

    if ~isempty(a1)
        accH1(k) = mean(a1, 'omitnan');
    else
        accH1(k) = NaN;
    end

    if ~isempty(a0)
        accH0(k) = mean(a0, 'omitnan');
    else
        accH0(k) = NaN;
    end
end

%% Print summary table
T = table(method, BER_mean, BER_std, H1_mean, H0_mean, gap_mean, accH1, accH0);
disp(T);

%% Figure 1: OTA trial-level comparison
figure('Color','w','Position',[100 80 950 850]);
tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% BER across trials
nexttile; hold on; grid on; box on;
for k = 1:nM
    plot(1:numel(allBER{k}), allBER{k}, '-o', 'LineWidth', 1.8, 'MarkerSize', 6);
end
ylabel('BER');
title('(a) OTA BER Across Trials');
legend(method, 'Location','best');
set(gca,'FontSize',12);

% H1/H0 scores
nexttile; hold on; grid on; box on;
x = 1:nM;
bar(x-0.16, H1_mean, 0.32);
bar(x+0.16, H0_mean, 0.32);
errorbar(x-0.16, H1_mean, H1_std, 'k.', 'LineWidth', 1.1);
errorbar(x+0.16, H0_mean, H0_std, 'k.', 'LineWidth', 1.1);
xticks(x); xticklabels(method);
ylabel('Mean authentication score');
title('(b) OTA H_1/H_0 Score Separation');
legend('H_1 legit','H_0 impostor','Location','best');
set(gca,'FontSize',12);

% Score gap
nexttile; hold on; grid on; box on;
bar(x, gap_mean);
errorbar(x, gap_mean, gap_std, 'k.', 'LineWidth', 1.2);
xticks(x); xticklabels(method);
ylabel('\Delta T = mean(T_{H1}) - mean(T_{H0})');
title('(c) OTA Authentication Score Gap');
set(gca,'FontSize',12);

sgtitle(tl,'OTA Authentication Comparison Across Methods','FontWeight','bold');

%% Figure 2: acceptance summary, if flags exist
figure('Color','w','Position',[120 120 780 460]);
hold on; grid on; box on;
bar(x-0.16, accH1, 0.32);
bar(x+0.16, accH0, 0.32);
xticks(x); xticklabels(method);
ylabel('Acceptance probability');
ylim([0 1]);
title('OTA Frame-Level Acceptance Summary');
legend('H_1 acceptance','H_0 acceptance','Location','best');
set(gca,'FontSize',12);

%% Save figures
saveas(figure(1), 'ota_compare_trials_scores.png');
saveas(figure(2), 'ota_acceptance_summary.png');

end

%% ------------------------------------------------------------------------
function v = get_field_vector(R, fieldName)

v = [];

if isstruct(R)
    if numel(R) > 1
        if isfield(R, fieldName)
            v = [R.(fieldName)];
        end
    else
        if isfield(R, fieldName)
            v = R.(fieldName);
        end
    end
end

v = v(:).';

end

%% ------------------------------------------------------------------------
function v = get_raw_metric_fallback(R, methodName, isH1)

v = [];

switch lower(string(methodName))
    case "permutation"
        if isH1
            rawField = 'Tp_H1';
        else
            rawField = 'Tp_H0';
        end

    case "im"
        if isH1
            rawField = 'Tim_H1';
        else
            rawField = 'Tim_H0';
        end

    case "dppa"
        if isH1
            rawField = 'Td_H1';
        else
            rawField = 'Td_H0';
        end

    otherwise
        return;
end

if numel(R) > 1
    vals = nan(1,numel(R));
    for i = 1:numel(R)
        if isfield(R(i), rawField)
            vals(i) = mean(R(i).(rawField), 'omitnan');
        end
    end
    v = vals;
else
    if isfield(R, rawField)
        raw = R.(rawField);
        if isvector(raw)
            v = mean(raw, 'omitnan');
        else
            v = mean(raw, 2, 'omitnan');
        end
    end
end

v = v(:).';

end