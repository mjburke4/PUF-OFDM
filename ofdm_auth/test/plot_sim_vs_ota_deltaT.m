function plot_sim_vs_ota_deltaT()

% Simulation vs OTA score-gap comparison
% DeltaT = mean(T_H1) - mean(T_H0)

clear; clc; close all;

%% -----------------------------
% 1) OTA values from latest runs
% ------------------------------
% Replace these with values from your saved S_perm/S_im/S_dppa if desired.

ota_perm = 0.04;   % from latest combined plot / permutation run, approx weak
ota_im   = 0.77;   % IM OTA score gap
ota_dppa = 0.74;   % DPPA OTA score gap

OTA = [ota_perm, ota_im, ota_dppa];

%% -----------------------------
% 2) Simulation values
% ------------------------------
% Fill these from your histogram or sim summary output:
% DeltaT_sim = mean(T_H1_sim) - mean(T_H0_sim)
load('results_sim.mat')

% sim_perm = 0.95;   % example placeholder
% sim_im   = 0.98;   % example placeholder
% sim_dppa = 0.45;   % example placeholder

% sim_perm = mean(results_sim.Tp_H1, 'omitnan')  - mean(results_sim.Tp_H0, 'omitnan');
% sim_im   = mean(results_sim.Tim_H1, 'omitnan') - mean(results_sim.Tim_H0, 'omitnan');
% sim_dppa = mean(results_sim.Td_H1, 'omitnan')  - mean(results_sim.Td_H0, 'omitnan');

sim_perm = results_sim.perm_deltaT;
sim_im   = results_sim.im_deltaT;
sim_dppa = results_sim.dppa_deltaT;

SIM = [sim_perm, sim_im, sim_dppa];

fprintf('Sim ΔT (Perm) = %.3f\n', sim_perm);
fprintf('Sim ΔT (IM)   = %.3f\n', sim_im);
fprintf('Sim ΔT (DPPA) = %.3f\n', sim_dppa);

%% -----------------------------
% 3) Plot
% ------------------------------
methods = categorical({'Permutation','IM','DPPA'});
methods = reordercats(methods, {'Permutation','IM','DPPA'});

Y = [SIM(:), OTA(:)];

figure('Color','w','Position',[120 120 760 460]);
bar(methods, Y, 'grouped');
grid on; box on;

ylabel('\DeltaT = mean(T_{H1}) - mean(T_{H0})');
title('Simulation vs OTA Authentication Score Separation');
legend('Simulation','OTA','Location','best');

ylim([0 1.1]);

set(gca,'FontSize',12,'LineWidth',1.0);

% Optional annotation
text(1, OTA(1)+0.05, 'receiver-limited', ...
    'HorizontalAlignment','center', 'FontSize',10);

exportgraphics(gcf, 'sim_vs_ota_deltaT.png', 'Resolution', 300);
savefig(gcf, 'sim_vs_ota_deltaT.fig');

end