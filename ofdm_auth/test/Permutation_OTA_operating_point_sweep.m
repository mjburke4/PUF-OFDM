%% ============================================
% Permutation OTA operating-point sweep
% Uses stored results(k).Tp_H1 and results(k).Tp_H0
%% ============================================

ok_idx = find([results.ok]);
if isempty(ok_idx)
    error('No successful trials found in results.');
end

tau_grid   = 0.0:0.05:3.0;     % permutation scores can exceed 1, so use wider range
alpha_grid = 0.40:0.05:0.9;
Lauth_grid = [8 16 32 64];

Ntau   = numel(tau_grid);
Nalpha = numel(alpha_grid);
NLauth = numel(Lauth_grid);

H1_acc_cube = zeros(Nalpha, Ntau, NLauth);
H0_acc_cube = zeros(Nalpha, Ntau, NLauth);
J_cube      = zeros(Nalpha, Ntau, NLauth);   % H1 - H0

for il = 1:NLauth
    Lauth = Lauth_grid(il);

    for ia = 1:Nalpha
        alpha = alpha_grid(ia);

        for it = 1:Ntau
            tau_p = tau_grid(it);

            accH1 = zeros(1, numel(ok_idx));
            accH0 = zeros(1, numel(ok_idx));

            for kk = 1:numel(ok_idx)
                r = results(ok_idx(kk));

                TpH1 = r.Tp_H1(:).';
                TpH0 = r.Tp_H0(:).';

                Lauth_eff = min([Lauth, numel(TpH1), numel(TpH0)]);

                TpH1a = TpH1(1:Lauth_eff);
                TpH0a = TpH0(1:Lauth_eff);

                VH1 = sum(TpH1a >= tau_p);
                VH0 = sum(TpH0a >= tau_p);

                reqVotes = ceil(alpha * Lauth_eff);

                accH1(kk) = (VH1 >= reqVotes);
                accH0(kk) = (VH0 >= reqVotes);
            end

            H1_acc_cube(ia,it,il) = mean(accH1);
            H0_acc_cube(ia,it,il) = mean(accH0);
            J_cube(ia,it,il)      = H1_acc_cube(ia,it,il) - H0_acc_cube(ia,it,il);
        end
    end
end

%% Find best operating point overall
bestJ = -inf;
best = struct();

for il = 1:NLauth
    [Jmax_this, idx] = max(J_cube(:,:,il), [], 'all', 'linear');
    if Jmax_this > bestJ
        bestJ = Jmax_this;
        [ia, it] = ind2sub([Nalpha, Ntau], idx);

        best.J = Jmax_this;
        best.alpha = alpha_grid(ia);
        best.tau_p = tau_grid(it);
        best.Lauth = Lauth_grid(il);
        best.H1acc = H1_acc_cube(ia,it,il);
        best.H0acc = H0_acc_cube(ia,it,il);
    end
end

fprintf('\n===== BEST PERM OTA OPERATING POINT =====\n');
fprintf('Best score J = H1-H0 = %.3f\n', best.J);
fprintf('tau_p  = %.3f\n', best.tau_p);
fprintf('alpha  = %.2f\n', best.alpha);
fprintf('Lauth  = %d\n', best.Lauth);
fprintf('H1 acceptance rate = %.3f\n', best.H1acc);
fprintf('H0 acceptance rate = %.3f\n', best.H0acc);

%% Plot heatmaps for each Lauth
for il = 1:NLauth
    figure;
    imagesc(tau_grid, alpha_grid, J_cube(:,:,il));
    set(gca,'YDir','normal');
    colorbar;
    xlabel('\tau_p');
    ylabel('\alpha');
    title(sprintf('Permutation OTA operating-point score: H1-H0 (Lauth=%d)', Lauth_grid(il)));
end

%% Also plot H1/H0 acceptance for the best Lauth
[~, il_best] = min(abs(Lauth_grid - best.Lauth));

figure;
imagesc(tau_grid, alpha_grid, H1_acc_cube(:,:,il_best));
set(gca,'YDir','normal');
colorbar;
xlabel('\tau_p');
ylabel('\alpha');
title(sprintf('Permutation OTA H1 acceptance (Lauth=%d)', Lauth_grid(il_best)));

figure;
imagesc(tau_grid, alpha_grid, H0_acc_cube(:,:,il_best));
set(gca,'YDir','normal');
colorbar;
xlabel('\tau_p');
ylabel('\alpha');
title(sprintf('Permutation OTA H0 acceptance (Lauth=%d)', Lauth_grid(il_best)));

ok_idx = find([results.ok]);
dvec = zeros(size(ok_idx));

for k = 1:numel(ok_idx)
    r = results(ok_idx(k));
    Lauth_eff = min(16, min(numel(r.Tp_H1), numel(r.Tp_H0)));
    dvec(k) = mean(r.Tp_H1(1:Lauth_eff)) - mean(r.Tp_H0(1:Lauth_eff));
end

figure;
plot(ok_idx, dvec, 'o-','LineWidth',1.5);
grid on;
xlabel('Trial');
ylabel('Mean(Tp\_H1 - Tp\_H0)');
title('Permutation OTA score gap by trial');