function smoke_test_multinode_Lsymbols()

cfg = cfg_default();

% Defaults
U = 10; m = 64; B = 64;
L = cfg.Lsym;

% Node population (must include: .A and distance fields)
nodes = init_node_population(U, m, 'seed', 7);

% Base challenge and per-frame nonce
C_seed = randi([0 1], 1, m);
nonce_frame  = randi([0 1], 1, 32);
nonce_replay = randi([0 1], 1, 32);   % stale/wrong nonce example (replay mismatch)

% Choose transmitting node
u_tx = 3;          % who actually transmits
u_expected = 3;    % who the slot expects (TDMA-style); set u_expected ~= u_tx to test spoofing

% Pathloss parameters (global)
%d0 = 1;            % meters
d0 = 10;            % meters
gamma = 2.7;       % pathloss exponent
SNRref_dB = 10;    % reference SNR at d0 (noise floor set from this)

% Thresholds (tune later)
tauP  = 0.05;
%tauP = 0.5 * median(Sp([1:u_tx-1,u_tx+1:end]) / L);

tauIM = 0.5;

% Accumulators per hypothesis user
Sp    = zeros(1,U);
Sim   = zeros(1,U);
Vperm = zeros(1,U);
Vim   = zeros(1,U);

cfg.chan.type = 'awgn';

% ---- Print node distances / gains (Step 4)
dist = arrayfun(@(n) n.distance_m, nodes);
gpow = (d0 ./ dist) .^ gamma;              % power gain
gdb  = 10*log10(gpow);

disp(table((1:U).', dist(:), gdb(:), 'VariableNames', {'User','distance_m','pathlossGain_dB'}));

% Fixed receiver noise floor based on reference SNR (Step 1 realism)
% We assume unit-average-power x_cp at amp=1. This is a good approximation.
noise_pow = 1 / (10^(SNRref_dB/10));

for ell = 1:L

    [X, Xp, meta] = gen_ofdm_symbol(cfg);
    % if ell == 1
    %     disp(Xp(meta.pilot_idx).');  % print pilot values
    % end

    % --- TX mapping (per-symbol)
    map_tx = verifier_mapping_apuf_sym(nodes(u_tx), C_seed, nonce_frame, ell, B, cfg.Npilots, cfg.K_active_pilots);
    [Xtilde, ~, ~] = apply_structure_operator_pilots(X, Xp, meta, map_tx.pilot_perm_local, map_tx.pilot_mask_local);

    % OFDM modulate
    x = ifft(Xtilde, cfg.Nfft);
    x_cp = [x(end-cfg.Ncp+1:end); x];

    % --- Reference signal power at d0 (amp = 1)
    Pref = mean(abs(x_cp).^2);

    % --- Set receiver noise floor so that SNRref holds at d0
    noise_pow = Pref / (10^(SNRref_dB/10));

    % Apply pathloss (amplitude scale)
    amp = sqrt((d0 / nodes(u_tx).distance_m)^gamma);

    xchan_in = amp * x_cp;

    % Channel
    [y_cp, ~] = channel_model(cfg, xchan_in);

     SNR_eff_dB = 10*log10((mean(abs(y_cp).^2)) / noise_pow);
    if ell == 1
        fprintf('Approx TX effective SNR at RX (ell=1): %.2f dB\n', SNR_eff_dB);
    end

    SNR_eff_pred = SNRref_dB + 10*log10((d0 / nodes(u_tx).distance_m)^gamma);
    SNR_eff_emp  = 10*log10(mean(abs(y_cp).^2) / noise_pow);
    
    if ell == 1
        fprintf('Pred TX SNR @ RX: %.2f dB | Empirical: %.2f dB\n', SNR_eff_pred, SNR_eff_emp);
    end


    % Add AWGN with fixed noise floor
    w = sqrt(noise_pow/2) * (randn(size(y_cp)) + 1j*randn(size(y_cp)));
    r_cp = y_cp + w;

    % FFT
    r = r_cp(cfg.Ncp+1:end);
    Y = fft(r, cfg.Nfft);

    % --- PERFECT CSI FOR PATHLOSS REALISM (do this first; LS later)
    % Channel model is AWGN => H = 1, but pathloss amp scales the whole link.
    % For general channel_model, you can have it return H (freq response); then multiply by amp.
    % Hhat = amp * ones(cfg.Nfft,1);   % since cfg.chan.type = 'awgn'
    % Z = Y ./ (Hhat + 1e-12);
    
    %Recommended for realism:
     %keep pathloss in the signal
    %fixed noise floor
    %equalize only the multipath channel (later), not the large-scale gain
    
    Z = Y;   % do not divide by amp

    % ---- Hypothesis loop (Metric A + B), with hypothesis-consistent LS channel estimation (Step 5)
    for u = 1:U
        % map_rx = verifier_mapping_apuf_sym(nodes(u), C_seed, nonce_frame, ell, B, cfg.Npilots, cfg.K_active_pilots);
        % [~, templ_rx, st_rx] = apply_structure_operator_pilots(X, Xp, meta, map_rx.pilot_perm_local, map_rx.pilot_mask_local);
        % 
        % % LS channel estimate consistent with predicted pilot structure
        % Hhat_u = estimate_channel_ls_pilots(Y, templ_rx, meta.pilot_idx, cfg.Nfft);
        % Z_u = Y ./ (Hhat_u + 1e-12);
        % 
        % % Metric A
        % tP = abs(metric_perm_correlator(Z_u, templ_rx));
        % Sp(u) = Sp(u) + tP;
        % Vperm(u) = Vperm(u) + (tP >= tauP);
        % 
        % % Metric B
        % tIM = metric_im_energysep_pilots(Z_u, meta.pilot_idx, st_rx.Atilde);
        % if isstruct(tIM), tIM = tIM.norm; end
        % Sim(u) = Sim(u) + tIM;
        % Vim(u) = Vim(u) + (tIM >= tauIM);
        map_rx = verifier_mapping_apuf_sym(nodes(u), C_seed, nonce_frame, ell, B, cfg.Npilots, cfg.K_active_pilots);
        [~, templ_rx, st_rx] = apply_structure_operator_pilots(X, Xp, meta, map_rx.pilot_perm_local, map_rx.pilot_mask_local);
        
        % Metric A on hypothesis template
        %tP = abs(metric_perm_correlator(Z, templ_rx));
        %tP = abs(metric_perm_correlator_pilots(Z, templ_rx, meta.pilot_idx));
        tP = abs(metric_perm_correlator_activepilots(Z, templ_rx, st_rx.Atilde));

        Sp(u) = Sp(u) + tP;
        Vperm(u) = Vperm(u) + (tP >= tauP);


        % Metric B on hypothesis active set
        tIM = metric_im_energysep_pilots(Z, meta.pilot_idx, st_rx.Atilde);
        if isstruct(tIM), tIM = tIM.norm; end
        Sim(u) = Sim(u) + tIM;
        Vim(u) = Vim(u) + (tIM >= tauIM);

    end

end

% Winners
[~, best_perm_soft] = max(Sp);
%[~, best_perm_vote] = max(Vperm);
if all(Vperm==0)
    fprintf('Perm vote: no hypothesis exceeded tauP\n');
else
    [~, best_perm_vote] = max(Vperm);
end

[~, best_im_soft]   = max(Sim);
[~, best_im_vote]   = max(Vim);

fprintf('Avg |Tperm| (TX): %.4f\n', Sp(u_tx)/L);
fprintf('Avg Tim (TX): %.4f\n', Sim(u_tx)/L);

fprintf('TX user u_tx=%d (expected=%d), L=%d, SNRref=%g dB at d0=%g m\n', u_tx, u_expected, L, SNRref_dB, d0);

T = table((1:U).', Sp.', Vperm.', Sim.', Vim.', ...
    'VariableNames', {'User','S_P_sum','V_P_votes','S_IM_sum','V_IM_votes'});
disp(T);

fprintf('Best (perm soft) user: %d\n', best_perm_soft);
fprintf('Best (perm vote) user: %d\n', best_perm_vote);
fprintf('Best (IM soft)   user: %d\n', best_im_soft);
fprintf('Best (IM vote)   user: %d\n', best_im_vote);

% ---- Quick “observe separation” summary (Step 5)
fprintf('TX distance=%.2f m, pathlossGain=%.2f dB\n', nodes(u_tx).distance_m, 10*log10((d0/nodes(u_tx).distance_m)^gamma));

% Show who is closest/farthest (sanity check)
[~, idx_close] = min(dist);
[~, idx_far]   = max(dist);
fprintf('Closest node: %d at %.2f m\n', idx_close, dist(idx_close));
fprintf('Farthest node: %d at %.2f m\n', idx_far, dist(idx_far));

% Show margin between best and second-best for each metric (separation behavior)
Sp_sorted = sort(Sp,'descend');
Sim_sorted = sort(Sim,'descend');
fprintf('Perm separation margin (best - 2nd): %.3f\n', Sp_sorted(1) - Sp_sorted(2));
fprintf('IM separation margin   (best - 2nd): %.3f\n', Sim_sorted(1) - Sim_sorted(2));

fprintf('Avg |Tperm| (TX user): %.4f\n', Sp(u_tx)/L);
fprintf('Avg Tim (TX user): %.4f\n', Sim(u_tx)/L);


end


function nonce2 = new_nonce(nbits)
    nonce2 = randi([0 1], 1, nbits);
end
