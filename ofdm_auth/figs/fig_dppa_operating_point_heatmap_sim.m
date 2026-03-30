function fig_dppa_operating_point_heatmap_sim()

cfg = cfg_default_dppa();
cfg.Lsym = 4;

U = 10;
nodes = init_node_population(U, cfg.m, 'seed', 7);
u0 = 3;

SNR_dB = 10;      % pick one representative SNR first
MC = 1000;        % can reduce to 300 first if slow

tau_grid   = 0.88:0.01:0.97;
alpha_grid = 0.60:0.05:0.90;

H1_acc_rate = zeros(numel(alpha_grid), numel(tau_grid));
H0_acc_rate = zeros(numel(alpha_grid), numel(tau_grid));
Jscore      = zeros(numel(alpha_grid), numel(tau_grid));

for ia = 1:numel(alpha_grid)
    alpha_test = alpha_grid(ia);

    for it = 1:numel(tau_grid)
        tau_test = tau_grid(it);

        cfg_tmp = cfg;
        cfg_tmp.alpha = alpha_test;
        cfg_tmp.tauD  = tau_test;

        accH1 = 0;
        accH0 = 0;

        for t = 1:MC
            % Legit trial
            o1 = run_one_trial_dppa(cfg_tmp, nodes, u0, SNR_dB, true);
            accH1 = accH1 + o1.accept_dppa;

            % Impostor trial
            o0 = run_one_trial_dppa(cfg_tmp, nodes, u0, SNR_dB, false);
            accH0 = accH0 + o0.accept_dppa;
        end

        H1_acc_rate(ia,it) = accH1 / MC;
        H0_acc_rate(ia,it) = accH0 / MC;
        Jscore(ia,it)      = H1_acc_rate(ia,it) - H0_acc_rate(ia,it);
    end
end

% Main score heatmap
figure;
imagesc(tau_grid, alpha_grid, Jscore);
set(gca,'YDir','normal');
colorbar;
xlabel('\tau_D');
ylabel('\alpha');
title(sprintf('DPPA sim operating-point score at SNR=%d dB: H1 acceptance - H0 acceptance', SNR_dB));

% Optional separate heatmaps
figure;
imagesc(tau_grid, alpha_grid, H1_acc_rate);
set(gca,'YDir','normal');
colorbar;
xlabel('\tau_D');
ylabel('\alpha');
title(sprintf('DPPA sim H1 acceptance rate at SNR=%d dB', SNR_dB));

figure;
imagesc(tau_grid, alpha_grid, H0_acc_rate);
set(gca,'YDir','normal');
colorbar;
xlabel('\tau_D');
ylabel('\alpha');
title(sprintf('DPPA sim H0 acceptance rate at SNR=%d dB', SNR_dB));

end