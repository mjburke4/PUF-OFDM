function out = test_one_symbol_smoke(SNRdB, is_legit)
% test_one_symbol_smoke  Quick end-to-end sanity check for OFDM auth metrics.
%   out = test_one_symbol_smoke(SNRdB, is_legit)
%
% Inputs:
%   SNRdB    : scalar
%   is_legit : true -> receiver uses correct mapping; false -> receiver uses wrong mapping
%
% Output struct:
%   out.Tperm_abs : |T^(P)| correlator magnitude
%   out.Tim       : T^(IM) energy-separation statistic
%   out.debug     : misc info (active set, pilots, etc.)

cfg = cfg_default();

% ---- 1) Generate logical symbol (data + pilots) and pilot-only vector
[X, Xp, meta] = gen_ofdm_symbol(cfg);

% ---- 2) Choose a "PUF-derived seed" (deterministic for repeatability)
% In real sim: seed_u32 = hash(PUF_response || nonce || symbol_index || user_id)
seed_legit = uint32(123456789);

% Transmit mapping (always legit at TX)
[Xtilde, templ_tx, st_tx] = apply_structure_operator(cfg, X, Xp, seed_legit, meta);

% ---- 3) OFDM modulation (IFFT + CP)
x = ifft(Xtilde, cfg.Nfft);
x_cp = [x(end-cfg.Ncp+1:end); x];

% ---- 4) Channel (AWGN baseline)
cfg.chan.type = 'awgn';
[y_cp, H] = channel_model(cfg, x_cp); % H is FFT of channel taps (size Nfft)

% Add AWGN in time domain
sig_pow = mean(abs(y_cp).^2);
SNRlin = 10^(SNRdB/10);
noise_pow = sig_pow / SNRlin;
w = sqrt(noise_pow/2) * (randn(size(y_cp)) + 1j*randn(size(y_cp)));
r_cp = y_cp + w;

% ---- 5) Receiver: remove CP + FFT
r = r_cp(cfg.Ncp+1:end);
Y = fft(r, cfg.Nfft);

% ---- 6) Perfect equalization (since channel is AWGN, H ~ 1 anyway)
Hhat = estimate_channel(cfg, H, [], [], meta.pilot_idx);
Z = Y ./ (Hhat + 1e-12);

% ---- 7) Receiver hypothesizes mapping (legit or not)
if is_legit
    seed_rx = seed_legit;
else
    seed_rx = uint32(987654321); % wrong mapping hypothesis
end
[~, templ_rx, st_rx] = apply_structure_operator(cfg, X, Xp, seed_rx, meta);

% ---- 8) Metric A: pilot-template correlator
Tperm = metric_perm_correlator(Z, templ_rx);

% ---- 9) Metric B: energy separation (using predicted active pilot set)
%Tim = metric_im_energysep(Z, st_rx.Atilde, cfg.Nfft);
Tim = metric_im_energysep_pilots(Z, meta.pilot_idx, st_rx.Atilde);

% Outputs
out.Tperm_abs = abs(Tperm);
out.Tim = Tim;
out.debug.st_tx = st_tx;
out.debug.st_rx = st_rx;
out.debug.meta  = meta;
out.debug.SNRdB = SNRdB;
out.debug.is_legit = is_legit;
end
