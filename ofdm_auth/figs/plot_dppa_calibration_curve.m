function plot_dppa_calibration_curve()

cfgD = cfg_default_dppa();
cfgD.Lsym  = 16;
cfgD.alpha = 0.75;
cfgD.adv_mode = 'impostor_othernode';

U = 10;
m = cfgD.m;
nodes = init_node_population(U, m, 'seed', 7);

u0 = 3;
SNR_cal = 15; %15;

MC = 1000;   % increase if noisy

tauD_grid = linspace(0.5, 0.999, 100);

FAR = zeros(size(tauD_grid));
FRR = zeros(size(tauD_grid));

for ti = 1:numel(tauD_grid)
    cfgD.tauD = tauD_grid(ti);

    acc0 = 0;
    rej1 = 0;

    for t = 1:MC
        % Impostor
        o0 = run_one_trial_dppa(cfgD, nodes, u0, SNR_cal, false);
        acc0 = acc0 + o0.accept_dppa;

        % Legit
        o1 = run_one_trial_dppa(cfgD, nodes, u0, SNR_cal, true);
        rej1 = rej1 + (~o1.accept_dppa);
    end

    FAR(ti) = acc0 / MC;
    FRR(ti) = rej1 / MC;
end

figure; hold on; grid on;
plot(tauD_grid, FAR, 'LineWidth', 2);
plot(tauD_grid, FRR, 'LineWidth', 2);
xlabel('\tau_D');
ylabel('Probability');
legend('FAR','FRR');
title('DPPA Calibration Curves at 15 dB');

set(gca,'YScale','log');

end