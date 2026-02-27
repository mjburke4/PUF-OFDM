function fig_far_frr_vs_tau()

cfg = cfg_default();
cfg.chan.type = 'awgn';

% Make curves visible
cfg.Lsym  = 8;
cfg.alpha = 0.75;

U = 10; m = 64;
nodes = init_node_population(U, m, 'seed', 7);

% Fix all nodes at reference distance so SNRref is true SNR
d0 = 10;
for u = 1:U
    nodes(u).distance_m = d0;
end

u0 = 3;

SNR_list = [-5 5 15];     % SNRs to compare
MC = 5000;

%% ================= Permutation Plot =====================

% tauP_grid = linspace(0, 1, 200);
% 
% figure; hold on; grid on;
% 
% for s = 1:length(SNR_list)
%     snr = SNR_list(s);
% 
%     FAR = zeros(size(tauP_grid));
%     FRR = zeros(size(tauP_grid));
% 
%     for ti = 1:length(tauP_grid)
%         cfg.tauP = tauP_grid(ti);
% 
%         % Legit trials -> FRR
%         rej = 0;
%         for t = 1:MC
%             o = run_one_trial_auth(cfg, nodes, u0, snr, true);
%             rej = rej + (~o.accept_perm);
%         end
%         FRR(ti) = rej/MC;
% 
%         % Impostor trials -> FAR
%         acc = 0;
%         for t = 1:MC
%             o = run_one_trial_auth(cfg, nodes, u0, snr, false);
%             acc = acc + (o.accept_perm);
%         end
%         FAR(ti) = acc/MC;
%     end
% 
%     FAR_s = movmean(FAR, 5);
%     FRR_s = movmean(FRR, 5);
%     plot(tauP_grid, FAR_s, '-', 'LineWidth', 3);
%     plot(tauP_grid, FRR_s, '--', 'LineWidth', 3);
% end
% 
% xlabel('\tau_P');
% ylabel('Probability');
% %title('Permutation Authentication: FAR/FRR vs \tau_P');
% title(sprintf('Permutation Authentication: FAR/FRR vs \\tau_P (L=%d, \\alpha=%.2f)', ...
%     cfg.Lsym, cfg.alpha));
% 
% legendStrings = {};
% set(gca,'YScale','log')
% for s = 1:length(SNR_list)
%     legendStrings{end+1} = sprintf('FAR (SNR= %d dB)', SNR_list(s));
%     legendStrings{end+1} = sprintf('FRR (SNR= %d dB)', SNR_list(s));
% end
% legend(legendStrings, 'Location','best');


%% ================= IM Plot =====================

%tauIM_grid = linspace(0, 0.4, 20);
tauIM_grid = linspace(0, 1, 200);

figure; hold on; grid on;

for s = 1:length(SNR_list)
    snr = SNR_list(s);

    FAR = zeros(size(tauIM_grid));
    FRR = zeros(size(tauIM_grid));

    for ti = 1:length(tauIM_grid)
        cfg.tauIM = tauIM_grid(ti);

        % Legit trials -> FRR
        rej = 0;
        for t = 1:MC
            o = run_one_trial_auth(cfg, nodes, u0, snr, true);
            rej = rej + (~o.accept_im);
        end
        FRR(ti) = rej/MC;

        % Impostor trials -> FAR
        acc = 0;
        for t = 1:MC
            o = run_one_trial_auth(cfg, nodes, u0, snr, false);
             acc = acc + (o.accept_im);
        end
        FAR(ti) = acc/MC;
    end
    
    FAR_s = movmean(FAR, 5);
    FRR_s = movmean(FRR, 5);
    plot(tauIM_grid, FAR_s, '-', 'LineWidth', 3);
    plot(tauIM_grid, FRR_s, '--', 'LineWidth', 3);
end

xlabel('\tau_{IM}');
ylabel('Probability');
%title('IM Authentication: FAR/FRR vs \tau_{IM}');
title(sprintf('IM Authentication: FAR/FRR vs \\tau_{IM} (L=%d, \\alpha=%.2f)', ...
    cfg.Lsym, cfg.alpha));

legendStrings = {};
set(gca,'YScale','log')
for s = 1:length(SNR_list)
    legendStrings{end+1} = sprintf('FAR (SNR= %d dB)', SNR_list(s));
    legendStrings{end+1} = sprintf('FRR (SNR= %d dB)', SNR_list(s));
end
legend(legendStrings, 'Location','best');

end