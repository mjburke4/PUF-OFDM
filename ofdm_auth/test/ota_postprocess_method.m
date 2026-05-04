function S = ota_postprocess_method(results, methodName, scoreH1Field, scoreH0Field, tau_grid, alpha_grid)

ok_idx = find([results.ok]);

if isempty(ok_idx)
    error('No successful trials.');
end

ber_vec = [results(ok_idx).ber];

H1mean = zeros(1,numel(ok_idx));
H0mean = zeros(1,numel(ok_idx));

allH1 = [];
allH0 = [];

for n = 1:numel(ok_idx)
    r = results(ok_idx(n));

    x1 = r.(scoreH1Field);
    x0 = r.(scoreH0Field);

    H1mean(n) = mean(x1,'omitnan');
    H0mean(n) = mean(x0,'omitnan');

    allH1 = [allH1, x1(:).'];
    allH0 = [allH0, x0(:).'];
end

% Sweep tau and alpha using frame-level voting
H1_acc_rate = zeros(numel(alpha_grid), numel(tau_grid));
H0_acc_rate = zeros(numel(alpha_grid), numel(tau_grid));
Jscore      = zeros(numel(alpha_grid), numel(tau_grid));

for ia = 1:numel(alpha_grid)
    alpha_test = alpha_grid(ia);

    for it = 1:numel(tau_grid)
        tau_test = tau_grid(it);

        accH1 = false(1,numel(ok_idx));
        accH0 = false(1,numel(ok_idx));

        for n = 1:numel(ok_idx)
            r = results(ok_idx(n));

            x1 = r.(scoreH1Field);
            x0 = r.(scoreH0Field);

            L1 = numel(x1);
            L0 = numel(x0);

            VH1 = sum(x1 >= tau_test);
            VH0 = sum(x0 >= tau_test);

            accH1(n) = VH1 >= ceil(alpha_test * L1);
            accH0(n) = VH0 >= ceil(alpha_test * L0);
        end

        H1_acc_rate(ia,it) = mean(accH1);
        H0_acc_rate(ia,it) = mean(accH0);
        Jscore(ia,it) = H1_acc_rate(ia,it) - H0_acc_rate(ia,it);
    end
end

% Choose best operating point:
% maximize H1 acceptance subject to H0 acceptance <= targetFAR
targetFAR_ota = 0.05;

valid = H0_acc_rate <= targetFAR_ota;

if any(valid(:))
    metric = H1_acc_rate;
    metric(~valid) = -inf;
    [~,idx] = max(metric(:));
else
    warning('%s: no operating point met target FAR. Using max J-score.', methodName);
    [~,idx] = max(Jscore(:));
end

[bestAlphaIdx,bestTauIdx] = ind2sub(size(Jscore), idx);

best_alpha = alpha_grid(bestAlphaIdx);
best_tau   = tau_grid(bestTauIdx);

best_H1acc = H1_acc_rate(bestAlphaIdx,bestTauIdx);
best_H0acc = H0_acc_rate(bestAlphaIdx,bestTauIdx);

fprintf('\n===== %s OTA SUMMARY OVER %d SUCCESSFUL TRIALS =====\n', methodName, numel(ok_idx));
fprintf('Mean BER              = %.3e\n', mean(ber_vec,'omitnan'));
fprintf('Median BER            = %.3e\n', median(ber_vec,'omitnan'));
fprintf('Mean T_H1             = %.3f\n', mean(H1mean,'omitnan'));
fprintf('Mean T_H0             = %.3f\n', mean(H0mean,'omitnan'));
fprintf('Mean score gap        = %.3f\n', mean(H1mean-H0mean,'omitnan'));
fprintf('Best tau              = %.4f\n', best_tau);
fprintf('Best alpha            = %.3f\n', best_alpha);
fprintf('H1 acceptance rate    = %.3f\n', best_H1acc);
fprintf('H0 acceptance rate    = %.3f\n', best_H0acc);

% Table
% T = table( ...
%     [results.trial].', ...
%     [results.fc].'/1e6, ...
%     [results.ok].', ...
%     [results.ber].', ...
%     H1mean_all(results, ok_idx, scoreH1Field).', ...
%     H0mean_all(results, ok_idx, scoreH0Field).', ...
%     'VariableNames', {'trial','fc_MHz','ok','BER','T_H1_mean','T_H0_mean'});

T = table( ...
    [results(ok_idx).trial].', ...
    [results(ok_idx).fc].'/1e6, ...
    [results(ok_idx).ok].', ...
    [results(ok_idx).ber].', ...
    H1mean.', ...
    H0mean.', ...
    'VariableNames', {'trial','fc_MHz','ok','BER','T_H1_mean','T_H0_mean'} );

disp(T);

disp(T);

% Figure 1: BER
figure('Color','w');
plot(ok_idx, ber_vec, 'o-','LineWidth',1.5);
grid on; box on;
xlabel('Trial index');
ylabel('BER');
title(sprintf('%s: BER across OTA trials', methodName));

% Figure 2: score separation
figure('Color','w');
plot(ok_idx, H1mean, 'o-','LineWidth',1.5); hold on;
plot(ok_idx, H0mean, 's-','LineWidth',1.5);
grid on; box on;
legend('H_1 legit','H_0 impostor','Location','best');
xlabel('Trial index');
ylabel('Mean authentication score');
title(sprintf('%s: H_1/H_0 score separation', methodName));

% Figure 3: per-symbol pass fraction vs tau
passFrac_H1 = zeros(size(tau_grid));
passFrac_H0 = zeros(size(tau_grid));

for k = 1:numel(tau_grid)
    passFrac_H1(k) = mean(allH1 >= tau_grid(k));
    passFrac_H0(k) = mean(allH0 >= tau_grid(k));
end

figure('Color','w');
plot(tau_grid, passFrac_H1, 'LineWidth',1.7); hold on;
plot(tau_grid, passFrac_H0, 'LineWidth',1.7);
grid on; box on;
xlabel('\tau');
ylabel('Per-symbol pass fraction');
legend('H_1','H_0','Location','best');
title(sprintf('%s: per-symbol pass fraction vs threshold', methodName));

% Figure 4: operating point heatmap
figure('Color','w');
imagesc(tau_grid, alpha_grid, Jscore);
axis xy; colorbar;
xlabel('\tau');
ylabel('\alpha');
title(sprintf('%s: J = H_1 acceptance - H_0 acceptance', methodName));
hold on;
plot(best_tau, best_alpha, 'rx', 'LineWidth',2, 'MarkerSize',12);

% Pack output
S.method = methodName;
S.ok_idx = ok_idx;
S.ber_vec = ber_vec;
S.H1mean = H1mean;
S.H0mean = H0mean;
S.gap = H1mean - H0mean;
S.tau_grid = tau_grid;
S.alpha_grid = alpha_grid;
S.H1_acc_rate = H1_acc_rate;
S.H0_acc_rate = H0_acc_rate;
S.Jscore = Jscore;
S.best_tau = best_tau;
S.best_alpha = best_alpha;
S.best_H1acc = best_H1acc;
S.best_H0acc = best_H0acc;

end

function v = H1mean_all(results, ok_idx, fieldName)
v = nan(1,numel(results));
for n = 1:numel(ok_idx)
    idx = ok_idx(n);
    v(idx) = mean(results(idx).(fieldName),'omitnan');
end
end

function v = H0mean_all(results, ok_idx, fieldName)
v = nan(1,numel(results));
for n = 1:numel(ok_idx)
    idx = ok_idx(n);
    v(idx) = mean(results(idx).(fieldName),'omitnan');
end
end