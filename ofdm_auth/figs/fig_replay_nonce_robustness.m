function fig_replay_nonce_robustness()

cfg = cfg_default();
cfg.chan.type = 'awgn';

% Choose operating point
cfg.Lsym  = 16;        % change to 32 later
cfg.alpha = 0.75;
cfg.tauP  = 0.180;    % use your tuned values
cfg.tauIM = 0.045;

cfg.tauIM = cfg.tauIM + 0.05;
%cfg.tauP  = 0.180 + 0.05;

accH = 0;
U = 10; m = 64;
nodes = init_node_population(U, m, 'seed', 7);

% Fix distance so SNRref is meaningful
d0 = 10;
for u = 1:U
    nodes(u).distance_m = d0;
end

u0 = 3;
SNR_list = [-5 5 15 25];     % two points is usually enough
MC = 200000;

% Outputs
legit_FRR_P = zeros(size(SNR_list));
legit_FRR_IM = zeros(size(SNR_list));
replay_FAR_P = zeros(size(SNR_list));
replay_FAR_IM = zeros(size(SNR_list));
replay_FAR_H = zeros(size(SNR_list));

for si = 1:numel(SNR_list)
    snr = SNR_list(si);

    % Legit: fresh nonce
    rejP = 0; rejIM = 0;
    % for t = 1:MC
    %     o = run_one_trial_auth(cfg, nodes, u0, snr, true, 0); % p_flip=0
    %     rejP  = rejP  + (~o.accept_perm);
    %     rejIM = rejIM + (~o.accept_im);
    % end
    % legit_FRR_P(si)  = rejP/MC;
    % legit_FRR_IM(si) = rejIM/MC;

    % Replay: attacker replays old frame while verifier expects new nonce.
    % We simulate this by generating a "TX mapping" with stale nonce but building
    % the verifier hypothesis with a fresh nonce.
    accP = 0; accIM = 0; accH = 0;

    for t = 1:MC
        % One frame: create two different nonces
        nonce_old = randi([0 1], 1, 32);
        nonce_new = randi([0 1], 1, 32);

        % Run one replay trial (helper below)
        out = run_one_trial_replay(cfg, nodes, u0, snr, nonce_old, nonce_new);

        accP  = accP  + out.accept_perm;
        accIM = accIM + out.accept_im;
        accH = accH + out.accept_hybrid;
    end

    replay_FAR_P(si)  = accP/MC;
    replay_FAR_IM(si) = accIM/MC;
    replay_FAR_H(si) = accH / MC;

    fprintf('SNR=%g dB: Legit FRR(P)=%.4f FRR(IM)=%.4f | Replay FAR(P)=%.4f FAR(IM)=%.4f FAR(Hybrid)=%.4f\n', ...
        snr, legit_FRR_P(si), legit_FRR_IM(si), replay_FAR_P(si), replay_FAR_IM(si), replay_FAR_H(si));
end

% Plot bar chart (compact)
% figure; grid on; hold on;
% x = 1:numel(SNR_list);
% 
% bar(x-0.15, legit_FRR_P, 0.25);
% bar(x+0.15, replay_FAR_P, 0.25);
% set(gca,'XTick',x,'XTickLabel',string(SNR_list));
% xlabel('SNR_{ref} (dB)');
% ylabel('Probability');
% title(sprintf('Permutation: Legit FRR vs Replay FAR (L=%d, \\alpha=%.2f)', cfg.Lsym, cfg.alpha));
% legend('Legit FRR','Replay FAR', 'Location','best');
% 
% figure; grid on; hold on;
% bar(x-0.15, legit_FRR_IM, 0.25);
% bar(x+0.15, replay_FAR_IM, 0.25);
% set(gca,'XTick',x,'XTickLabel',string(SNR_list));
% xlabel('SNR_{ref} (dB)');
% ylabel('Probability');
% title(sprintf('IM: Legit FRR vs Replay FAR (L=%d, \\alpha=%.2f)', cfg.Lsym, cfg.alpha));
% legend('Legit FRR','Replay FAR', 'Location','best');

% --- Combined plot: Legit FRR vs Replay FAR for Perm and IM on one figure
figure; hold on; grid on;

x = 1:numel(SNR_list);

% Plot replay FAR curves
plot(x, max(replay_FAR_P, 1/MC), '-o', 'LineWidth', 2);
plot(x, max(replay_FAR_IM, 1/MC), '-s', 'LineWidth', 2);
plot(x, max(replay_FAR_H, 1/MC), '-d', 'LineWidth', 2);

% Optional: plot legit FRR as faint dashed lines
%plot(x, max(legit_FRR_P, 1/MC), '--o', 'LineWidth', 1);
%plot(x, max(legit_FRR_IM, 1/MC), '--s', 'LineWidth', 1);

set(gca,'XTick',x,'XTickLabel',string(SNR_list));
set(gca,'YScale','log');
ylim([1e-6 1]);

xlabel('SNR_{ref} (dB)');
ylabel('Probability (log scale)');
title(sprintf('Replay Robustness: Perm vs IM vs Hybrid (L=%d, \\alpha=%.2f)', ...
    cfg.Lsym, cfg.alpha));
% 
% legend('Replay FAR (Perm)', ...
%        'Replay FAR (IM)', ...
%        'Replay FAR (Hybrid)', ...
%        'Legit FRR (Perm)', ...
%        'Legit FRR (IM)', ...
%        'Location','best');

legend('Replay FAR (Perm)', ...
       'Replay FAR (IM)', ...
       'Replay FAR (Hybrid)', ...
       'Location','best');

% figure; hold on; grid on;
% 
% x = 1:numel(SNR_list);
% 
% % Bars: Replay FAR (Perm vs IM)
% bar(x-0.20, replay_FAR_P, 0.18);
% bar(x+0.00, replay_FAR_IM, 0.18);
% 
% % Markers: Legit FRR (Perm vs IM) as points (usually near zero)
% %plot(x-0.20, legit_FRR_P, 'ks', 'MarkerFaceColor','k');
% %plot(x+0.00, legit_FRR_IM, 'kd', 'MarkerFaceColor','k');
% 
% set(gca,'XTick',x,'XTickLabel',string(SNR_list));
% set(gca,'YScale','log');
% ylim([1e-4 1]);
% xlabel('SNR_{ref} (dB)');
% ylabel('Probability');
% title(sprintf('Replay Robustness: Perm vs IM (L=%d, \\alpha=%.2f)', cfg.Lsym, cfg.alpha));
% 
% legend('Replay FAR (Perm)','Replay FAR (IM)', 'Legit FRR (Perm)','Legit FRR (IM)', ...
%     'Location','best');
end


% ==========================================================
function out = run_one_trial_replay(cfg, nodes, u0, SNRref_dB, nonce_old, nonce_new)
% Replay trial:
% - TX uses stale nonce_old to build structure
% - Verifier expects nonce_new to build hypothesis/template

L = cfg.Lsym;
m = nodes(u0).m;
B = 64;

% Use same base challenge for both; only nonce differs
C_seed = randi([0 1], 1, m);

Vperm = 0;
Vim   = 0;

for ell = 1:L
    [X, Xp, meta] = gen_ofdm_symbol(cfg);

    % ---- TX uses OLD nonce
    map_tx = verifier_mapping_apuf_sym(nodes(u0), C_seed, nonce_old, ell, B, cfg.Npilots, cfg.K_active_pilots);
    [Xtilde, ~, ~] = apply_structure_operator_pilots(X, Xp, meta, map_tx.pilot_perm_local, map_tx.pilot_mask_local);

    x = ifft(Xtilde, cfg.Nfft);
    x_cp = [x(end-cfg.Ncp+1:end); x];

    Pref = mean(abs(x_cp).^2);
    noise_pow = Pref / (10^(SNRref_dB/10));

    [y_cp, ~] = channel_model(cfg, x_cp);
    w = sqrt(noise_pow/2) * (randn(size(y_cp)) + 1j*randn(size(y_cp)));
    r_cp = y_cp + w;

    r = r_cp(cfg.Ncp+1:end);
    Y = fft(r, cfg.Nfft);
    Z = Y;

    % ---- Verifier expects NEW nonce
    map_rx = verifier_mapping_apuf_sym(nodes(u0), C_seed, nonce_new, ell, B, cfg.Npilots, cfg.K_active_pilots);
    [~, templ_rx, st_rx] = apply_structure_operator_pilots(X, Xp, meta, map_rx.pilot_perm_local, map_rx.pilot_mask_local);

    % Perm metric (use your current choice)
    Tc = metric_perm_correlator_activepilots(Z, templ_rx, st_rx.Atilde);
    tP = real(Tc);

    % IM metric
    tIM = metric_im_energysep_pilots(Z, meta.pilot_idx, st_rx.Atilde);

    Vperm = Vperm + (tP >= cfg.tauP);
    Vim   = Vim   + (tIM >= cfg.tauIM);
end

need = ceil(cfg.alpha * cfg.Lsym);

out.Vperm = Vperm;
out.Vim   = Vim;
out.accept_perm = (Vperm >= need);
out.accept_im   = (Vim   >= need);
out.accept_hybrid = out.accept_perm && out.accept_im;
end