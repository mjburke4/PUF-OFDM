function fig_hist_metrics_H1_H0()

cfg = cfg_default();
cfg.chan.type = 'awgn';
cfg.Lsym  = 8;
cfg.alpha = 0.75;

% Keep it per-symbol: L not used for the histogram, but we’ll reuse your trial function style
U = 10; m = 64;
nodes = init_node_population(U, m, 'seed', 7);

% Fix distance so SNRref is actual receive SNR at d0
d0 = 10;
for u = 1:U
    nodes(u).distance_m = d0;
end

u0 = 3;

% Choose SNR points to illustrate (one or two)
SNR_list = [5];     % change to just [10] if you want one figure

Nsamp = 200000;         % per hypothesis, per SNR (increase for smoother hist)
B = 64;

% For reproducibility
rng(1);

for si = 1:numel(SNR_list)
    SNRref_dB = SNR_list(si);

    % Storage
    Tp_H1  = zeros(1, Nsamp);
    Tp_H0  = zeros(1, Nsamp);
    Tim_H1 = zeros(1, Nsamp);
    Tim_H0 = zeros(1, Nsamp);

    % Shared challenge/nonce per sample? In reality it changes per frame.
    % For histogram, we can refresh per sample to avoid accidental structure reuse.
    for n = 1:Nsamp
        % One OFDM symbol index (ell) – include it so mapping evolves
        ell = 1;  % doesn’t matter if nonce+ell mixing is used; keep fixed for per-symbol view

        % Fresh frame challenge/nonce each sample
        C_seed = randi([0 1], 1, m);
        nonce_frame = randi([0 1], 1, 32);

        % Generate logical symbol
        [X, Xp, meta] = gen_ofdm_symbol(cfg);

        % ---------------- H1: Legit TX uses u0 mapping ----------------
        map_tx = verifier_mapping_apuf_sym(nodes(u0), C_seed, nonce_frame, ell, B, cfg.Npilots, cfg.K_active_pilots);
        [Xtilde, ~, ~] = apply_structure_operator_pilots(X, Xp, meta, map_tx.pilot_perm_local, map_tx.pilot_mask_local);

        % OFDM modulate
        x = ifft(Xtilde, cfg.Nfft);
        x_cp = [x(end-cfg.Ncp+1:end); x];

        % Noise calibration at d0
        Pref = mean(abs(x_cp).^2);
        noise_pow = Pref / (10^(SNRref_dB/10));

        % Channel
        [y_cp, ~] = channel_model(cfg, x_cp);
        w = sqrt(noise_pow/2) * (randn(size(y_cp)) + 1j*randn(size(y_cp)));
        r_cp = y_cp + w;

        % FFT
        r = r_cp(cfg.Ncp+1:end);
        Y = fft(r, cfg.Nfft);
        Z = Y;  % AWGN only

        % Verifier hypothesis template (u0)
        map_rx = verifier_mapping_apuf_sym(nodes(u0), C_seed, nonce_frame, ell, B, cfg.Npilots, cfg.K_active_pilots);
        [~, templ_rx, st_rx] = apply_structure_operator_pilots(X, Xp, meta, map_rx.pilot_perm_local, map_rx.pilot_mask_local);

        % Perm metric: use your coherence-sensitive statistic
        Tc = metric_perm_correlator_activepilots(Z, templ_rx, st_rx.Atilde);
        Tp_H1(n) = real(Tc);

        % IM metric
        Tim_H1(n) = metric_im_energysep_pilots(Z, meta.pilot_idx, st_rx.Atilde);

        % ---------------- H0: Impostor TX uses random mapping ----------------
        seed_u32 = uint32(randi([0, 2^32-1]));
        pilot_perm_local = perm_from_seed(seed_u32, cfg.Npilots);
        pilot_mask_local = mask_from_seed(seed_u32, cfg.Npilots, cfg.K_active_pilots);

        [Xtilde0, ~, ~] = apply_structure_operator_pilots(X, Xp, meta, pilot_perm_local, pilot_mask_local);

        x0 = ifft(Xtilde0, cfg.Nfft);
        x0_cp = [x0(end-cfg.Ncp+1:end); x0];

        Pref0 = mean(abs(x0_cp).^2);
        noise_pow0 = Pref0 / (10^(SNRref_dB/10));

        [y0_cp, ~] = channel_model(cfg, x0_cp);
        w0 = sqrt(noise_pow0/2) * (randn(size(y0_cp)) + 1j*randn(size(y0_cp)));
        r0_cp = y0_cp + w0;

        r0 = r0_cp(cfg.Ncp+1:end);
        Y0 = fft(r0, cfg.Nfft);
        Z0 = Y0;

        % Verifier still tests u0 hypothesis (templ_rx, st_rx already computed)
        Tc0 = metric_perm_correlator_activepilots(Z0, templ_rx, st_rx.Atilde);
        Tp_H0(n) = real(Tc0);

        Tim_H0(n) = metric_im_energysep_pilots(Z0, meta.pilot_idx, st_rx.Atilde);
    end

    % -------- Plot histograms --------
    figure; hold on; grid on;
    histogram(Tp_H0, 40, 'Normalization','probability');
    histogram(Tp_H1, 15, 'Normalization','probability');
    xlabel('T^{(P)} per-symbol statistic (real part)');
    ylabel('Probability');
    title(sprintf('Permutation Metric Histogram (H1 vs H0) at SNRref=%g dB (L=%d, \\alpha=%.2f)', SNRref_dB, cfg.Lsym, cfg.alpha));
    legend('H0 (impostor)','H1 (legit)', 'Location','best');

    figure; hold on; grid on;
    histogram(Tim_H0, 40, 'Normalization','probability');
    histogram(Tim_H1, 15, 'Normalization','probability');
    xlabel('T^{(IM)} per-symbol statistic');
    ylabel('Probability');
    title(sprintf('IM Metric Histogram (H1 vs H0) at SNRref=%g dB (L=%d, \\alpha=%.2f)', SNRref_dB, cfg.Lsym, cfg.alpha));
    legend('H0 (impostor)','H1 (legit)', 'Location','best');

    % Optional quick summary print
    fprintf('SNR=%g dB: Perm mean H1=%.3f H0=%.3f | IM mean H1=%.3f H0=%.3f\n', ...
        SNRref_dB, mean(Tp_H1), mean(Tp_H0), mean(Tim_H1), mean(Tim_H0));
end

end