function fig_far_frr_vs_tau_dppa()
addpath 'C:\Users\michael.j.burke10\OneDrive - US Army\Desktop\Grad school stuff\PUF based dynamic spreading code generation\simulations\PUF_Modulation_ofdm_Sim\PUF-OFDM\ofdm_auth\lib'
addpath 'C:\Users\michael.j.burke10\OneDrive - US Army\Desktop\Grad school stuff\PUF based dynamic spreading code generation\simulations\PUF_Modulation_ofdm_Sim\PUF-OFDM\ofdm_auth\cfg'
addpath 'C:\Users\michael.j.burke10\OneDrive - US Army\Desktop\Grad school stuff\PUF based dynamic spreading code generation\simulations\PUF_Modulation_ofdm_Sim\PUF-OFDM\ofdm_auth\sim'
cfg = cfg_default_dppa();

cfg.Lsym = 16;

cfg.alpha = 0.75;

U = 10;
nodes = init_node_population(U, cfg.m, 'seed', 7);
u0 = 3;

SNR_list = [-5 5 15];
MC = 6000;

tauD_grid = linspace(0.5, 1, 120);

figure; hold on; grid on;

for s = 1:length(SNR_list)
    snr = SNR_list(s);

    FAR = zeros(size(tauD_grid));
    FRR = zeros(size(tauD_grid));

    for ti = 1:length(tauD_grid)
        cfg.tauD = tauD_grid(ti);

        rej = 0;
        for t = 1:MC
            o = run_one_trial_dppa(cfg, nodes, u0, snr, true);
            rej = rej + (~o.accept_dppa);
        end
        FRR(ti) = rej / MC;

        acc = 0;
        for t = 1:MC
            o = run_one_trial_dppa(cfg, nodes, u0, snr, false);
            acc = acc + o.accept_dppa;
        end
        FAR(ti) = acc / MC;
    end

    FAR_s = movmean(FAR, 5);
    FRR_s = movmean(FRR, 5);

    plot(tauD_grid, FAR_s, '-', 'LineWidth', 2.5);
    plot(tauD_grid, FRR_s, '--', 'LineWidth', 2.5);
end

xlabel('\tau_D');
ylabel('Probability');
title(sprintf('DPPA Authentication: FAR/FRR vs \\tau_D (L=%d, \\alpha=%.2f)', ...
    cfg.Lsym, cfg.alpha));
set(gca,'YScale','log');

legendStrings = {};
for s = 1:length(SNR_list)
    legendStrings{end+1} = sprintf('FAR (SNR=%d dB)', SNR_list(s));
    legendStrings{end+1} = sprintf('FRR (SNR=%d dB)', SNR_list(s));
end
legend(legendStrings, 'Location', 'best');

!git add .end