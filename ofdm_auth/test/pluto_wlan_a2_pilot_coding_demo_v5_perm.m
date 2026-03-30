clear; clc; close all;

Ntrials = 5;
fc_list = [350e6];
tau_p   = 0.35;
alpha   = 0.75;

fprintf('\n===== STARTING PERM OTA RUN =====\n');

for k = 1:Ntrials
    fc_this = fc_list(mod(k-1, numel(fc_list)) + 1);
    fprintf('\nTrial %d / %d  fc = %.3f MHz\n', k, Ntrials, fc_this/1e6);

    try
        res = run_one_perm_trial(fc_this, tau_p, alpha);
        fprintf('BER=%.3e, Tp_H1=%.3f, Tp_H0=%.3f, acc=%d/%d\n', ...
            res.ber, res.TpH1_mean, res.TpH0_mean, res.acceptH1, res.acceptH0);
    catch ME
        fprintf(2, 'Trial %d FAILED: %s\n', k, ME.message);
    end
end

%% ------------ Multi-trial settings ------------
Ntrials = 20;
fc_list = [350e6];
tau_p   = 0.35;    % starter value, tune later
alpha   = 0.75;

% Keep your current good operating point
phi_choices_deg = [-12, +12];

results = repmat(struct( ...
    'trial', [], ...
    'fc', [], ...
    'ber', [], ...
    'numErr', [], ...
    'Lbits', [], ...
    'TpH1_mean', [], ...
    'TpH0_mean', [], ...
    'VH1', [], ...
    'VH0', [], ...
    'Lsym', [], ...
    'acceptH1', [], ...
    'acceptH0', [], ...
    'pktStart', [], ...
    'coarseCFO', [], ...
    'fineCFO', [], ...
    'rxMaxAmp', [], ...
    'ok', false, ...
    'Tp_H1', [], ...
    'Tp_H0', [], ...
    'errmsg', ''), Ntrials, 1);
fprintf('\n===== STARTING DPPA MULTI-TRIAL RUN =====\n');

for k = 1:Ntrials
    fc_this = fc_list(mod(k-1, numel(fc_list)) + 1);

    fprintf('\n----------------------------------------\n');
    fprintf('Trial %d / %d   fc = %.3f MHz\n', k, Ntrials, fc_this/1e6);
    fprintf('----------------------------------------\n');

    try
    %res = run_one_dppa_trial(fc_this, tau_d, alpha, phi_choices_deg);
    res = run_one_perm_trial(fc_this, tau_p, alpha);
    
    % Fill preallocated struct field-by-field so MATLAB doesn't complain
    results(k).trial      = k;
    results(k).fc         = res.fc;
    results(k).ber        = res.ber;
    results(k).numErr     = res.numErr;
    results(k).Lbits      = res.Lbits;
    results(k).TpH1_mean  = res.TpH1_mean;
    results(k).TpH0_mean  = res.TpH0_mean;
    results(k).VH1        = res.VH1;
    results(k).VH0        = res.VH0;
    results(k).Lsym       = res.Lsym;
    results(k).acceptH1   = res.acceptH1;
    results(k).acceptH0   = res.acceptH0;
    results(k).pktStart   = res.pktStart;
    results(k).coarseCFO  = res.coarseCFO;
    results(k).fineCFO    = res.fineCFO;
    results(k).rxMaxAmp   = res.rxMaxAmp;
    results(k).ok         = true;
    results(k).Tp_H1      = res.Tp_H1;
    results(k).Tp_H0      = res.Tp_H0;
    results(k).errmsg     = '';

    fprintf('Trial %d summary: BER=%.3e, Tp_H1=%.3f, Tp_H0=%.3f, acc(H1/H0)=%d/%d\n', ...
    k, res.ber, res.TpH1_mean, res.TpH0_mean, res.acceptH1, res.acceptH0);

    catch ME
        results(k).trial = k;
        results(k).fc = fc_this;
        results(k).ok = false;
        results(k).errmsg = ME.message;

        fprintf(2, 'Trial %d FAILED: %s\n', k, ME.message);
    end
end

fprintf('\n===== MULTI-TRIAL RUN COMPLETE =====\n');

%% ------------ Summarize successful trials ------------
ok_idx = find([results.ok]);

if isempty(ok_idx)
    error('No successful trials completed.');
end

ber_vec      = [results(ok_idx).ber];
%TdH1_vec     = [results(ok_idx).TdH1_mean];
%TdH0_vec     = [results(ok_idx).TdH0_mean];
TpH1_vec = [results(ok_idx).TpH1_mean];
TpH0_vec = [results(ok_idx).TpH0_mean];
accH1_vec    = [results(ok_idx).acceptH1];
accH0_vec    = [results(ok_idx).acceptH0];
fc_ok        = [results(ok_idx).fc];

fprintf('\n===== SUMMARY OVER %d SUCCESSFUL TRIALS =====\n', numel(ok_idx));
fprintf('Mean BER              = %.3e\n', mean(ber_vec));
fprintf('Median BER            = %.3e\n', median(ber_vec));
fprintf('Mean Tp_H1            = %.3f\n', mean(TpH1_vec));
fprintf('Mean Tp_H0            = %.3f\n', mean(TpH0_vec));
fprintf('H1 acceptance rate    = %.3f\n', mean(accH1_vec));
fprintf('H0 acceptance rate    = %.3f\n', mean(accH0_vec));

%% ------------ Quick plots ------------
figure;
plot(ok_idx, ber_vec, 'o-','LineWidth',1.5);
grid on;
title('BER across successful DPPA trials');
xlabel('Trial index');
ylabel('BER');

figure;
plot(ok_idx, TdH1_vec, 'o-','LineWidth',1.5); hold on;
plot(ok_idx, TdH0_vec, 's-','LineWidth',1.5);
grid on;
legend('Mean T_{d,H1}','Mean T_{d,H0}','Location','best');
title('DPPA score separation across trials');
xlabel('Trial index');
ylabel('Mean differential score');

figure;
plot(ok_idx, dH1_vec, 'o-','LineWidth',1.5); hold on;
plot(ok_idx, dH0_vec, 's-','LineWidth',1.5);
grid on;
legend('Mean |dErr| H1','Mean |dErr| H0','Location','best');
title('Differential angle error across trials');
xlabel('Trial index');
ylabel('Degrees');

figure;
%gscatter(ber_vec, TdH1_vec - TdH0_vec, fc_ok/1e6);
fc_labels = categorical(compose('%.0f MHz', fc_ok/1e6));
gscatter(ber_vec, TpH1_vec - TpH0_vec, fc_labels);
grid on;
title('BER vs DPPA separation by center frequency');
xlabel('BER');
ylabel('\DeltaT_d = mean(T_{d,H1}) - mean(T_{d,H0})');
legend('Location','best');

%% ------------ Optional table for easy inspection ------------
T = table( ...
    [results.trial].', ...
    [results.fc].'/1e6, ...
    [results.ok].', ...
    [results.ber].', ...
    [results.TdH1_mean].', ...
    [results.TdH0_mean].', ...
    [results.dErrH1_mean_deg].', ...
    [results.dErrH0_mean_deg].', ...
    [results.acceptH1].', ...
    [results.acceptH0].', ...
    'VariableNames', {'trial','fc_MHz','ok','BER','TdH1_mean','TdH0_mean','dErrH1_deg','dErrH0_deg','acceptH1','acceptH0'});

disp(T);

%% ------------ DPPA operating-point sweep ------------
tau_grid   = 0.88:0.01:0.97;
alpha_grid = 0.60:0.05:0.90;

H1_acc_rate = zeros(numel(alpha_grid), numel(tau_grid));
H0_acc_rate = zeros(numel(alpha_grid), numel(tau_grid));
Jscore      = zeros(numel(alpha_grid), numel(tau_grid));   % simple separation score

ok_idx = find([results.ok]);

for ia = 1:numel(alpha_grid)
    alpha_test = alpha_grid(ia);

    for it = 1:numel(tau_grid)
        tau_test = tau_grid(it);

        accH1_tmp = zeros(numel(ok_idx),1);
        accH0_tmp = zeros(numel(ok_idx),1);

        for n = 1:numel(ok_idx)
            r = results(ok_idx(n));

            Ld = numel(r.Td_H1);
            VH1_test = sum(r.Td_H1 >= tau_test);
            VH0_test = sum(r.Td_H0 >= tau_test);

            accH1_tmp(n) = (VH1_test >= ceil(alpha_test * Ld));
            accH0_tmp(n) = (VH0_test >= ceil(alpha_test * Ld));
        end

        H1_acc_rate(ia,it) = mean(accH1_tmp);
        H0_acc_rate(ia,it) = mean(accH0_tmp);
        Jscore(ia,it)      = H1_acc_rate(ia,it) - H0_acc_rate(ia,it);
    end
end

figure;
plot(ok_idx, ber_vec, 'o-','LineWidth',1.5);
grid on;
xlabel('Trial index');
ylabel('BER');
title('BER across DPPA OTA trials at 940 MHz');

figure;
plot(ok_idx, TdH1_vec, 'o-','LineWidth',1.5); hold on;
plot(ok_idx, TdH0_vec, 's-','LineWidth',1.5);
grid on;
legend('Mean T_{d,H1}','Mean T_{d,H0}','Location','best');
xlabel('Trial index');
ylabel('Mean differential score');
title('DPPA score separation across OTA trials');

figure;
plot(ok_idx, dH1_vec, 'o-','LineWidth',1.5); hold on;
plot(ok_idx, dH0_vec, 's-','LineWidth',1.5);
grid on;
legend('Mean |dErr| H1','Mean |dErr| H0','Location','best');
xlabel('Trial index');
ylabel('Degrees');
title('Differential angle error across OTA trials');
tau_vec = 0.70:0.01:0.99;
passFrac_H1 = zeros(size(tau_vec));
passFrac_H0 = zeros(size(tau_vec));

allTdH1 = [];
allTdH0 = [];

for n = 1:numel(ok_idx)
    allTdH1 = [allTdH1, results(ok_idx(n)).Td_H1];
    allTdH0 = [allTdH0, results(ok_idx(n)).Td_H0];
end

for k = 1:numel(tau_vec)
    passFrac_H1(k) = mean(allTdH1 >= tau_vec(k));
    passFrac_H0(k) = mean(allTdH0 >= tau_vec(k));
end

figure;
plot(tau_vec, passFrac_H1, 'LineWidth',1.5); hold on;
plot(tau_vec, passFrac_H0, 'LineWidth',1.5);
grid on;
legend('H1 pass fraction','H0 pass fraction','Location','best');
xlabel('\tau_d');
ylabel('Fraction passing');
title('DPPA symbol-pass fraction vs threshold');

figure;
imagesc(tau_grid, alpha_grid, Jscore);
set(gca,'YDir','normal');
colorbar;
xlabel('\tau_d');
ylabel('\alpha');
title('DPPA operating-point score: H1 acceptance - H0 acceptance');

figure;
imagesc(tau_grid, alpha_grid, H1_acc_rate);
set(gca,'YDir','normal');
colorbar;
xlabel('\tau_d');
ylabel('\alpha');
title('H1 acceptance rate');

figure;
imagesc(tau_grid, alpha_grid, H0_acc_rate);
set(gca,'YDir','normal');
colorbar;
xlabel('\tau_d');
ylabel('\alpha');
title('H0 acceptance rate');

function res = run_one_dppa_trial(fc, tau_d, alpha, phi_choices_deg)

    %% ------------ USER SETTINGS ------------
    txID = 'usb:0';
    rxID = 'usb:1';

    fs = 20e6;
    txGain = -3;
    rxGain = 45;
    enablePilotCoding = true;

    samplesPerFrame = 300000;
    pauseTX = 0.2;

    %% ------------ WLAN config ------------
    cfgNonHT = wlanNonHTConfig;
    cfgNonHT.ChannelBandwidth = 'CBW20';
    cfgNonHT.MCS = 1;
    cfgNonHT.NumTransmitAntennas = 1;

    psduLen = 400;
    cfgNonHT.PSDULength = psduLen;

    psdu = randi([0 1], 8*cfgNonHT.PSDULength, 1, 'int8');

    ind = wlanFieldIndices(cfgNonHT);
    pktLen = ind.NonHTData(2);
    idleTime = 200e-6;

    %% ------------ PUF node setup ------------
    U = 1;
    u0 = 1;
    m = 64;
    B = 64;

    nodes = init_node_population(U, m, 'puf_bits', 128, 'seed', 7);
    node = nodes(u0);
    node_claimed = node;

    C_seed = randi([0 1], 1, node.m);
    nonce_frame = randi([0 1], 1, 32);

    %% ------------ Waveform generation ------------
    waveform_uncoded = wlanWaveformGenerator(psdu, cfgNonHT, 'IdleTime', idleTime);
    txPilotsBase = extractPayloadPilotsFromWaveform(waveform_uncoded, cfgNonHT);

    waveform = waveform_uncoded;
    if enablePilotCoding
        waveform = applyPilotSignCodingNonHT(waveform, cfgNonHT, node, C_seed, nonce_frame, B, phi_choices_deg);
    end

    guard = zeros(round(500e-6*fs),1);
    waveform = [guard; waveform; guard];
    waveform = 0.6 * waveform / max(abs(waveform)+1e-12);

    %% ------------ Pluto TX/RX objects ------------
    tx = sdrtx('Pluto', ...
        'RadioID', txID, ...
        'CenterFrequency', fc, ...
        'BasebandSampleRate', fs, ...
        'Gain', txGain);

    rx = sdrrx('Pluto', ...
        'RadioID', rxID, ...
        'CenterFrequency', fc, ...
        'BasebandSampleRate', fs, ...
        'SamplesPerFrame', samplesPerFrame, ...
        'GainSource', 'Manual', ...
        'Gain', rxGain, ...
        'OutputDataType', 'double');

    cleaner = onCleanup(@() cleanup_radios(tx, rx));

    %% ------------ TX/RX capture ------------
    transmitRepeat(tx, waveform);
    pause(pauseTX);

    maxTries = 8;
    pktStart = [];
    coarseCFO = NaN;
    y = [];

    for attempt = 1:maxTries
        y = rx();
        rxMaxAmp = max(abs(y));

        pktStartCand = wlanPacketDetect(y, cfgNonHT.ChannelBandwidth);
        if isempty(pktStartCand)
            fprintf('Attempt %d: no packet detected\n', attempt);
            continue;
        end

        pktStart = pktStartCand(1);
        yPktTry = y(pktStart:end);
        coarseCFOTry = wlanCoarseCFOEstimate(yPktTry, cfgNonHT.ChannelBandwidth);

        if abs(coarseCFOTry) > 50e3
            fprintf('Attempt %d: bad packet detect (CFO %.1f Hz)\n', attempt, coarseCFOTry);
            continue;
        end

        coarseCFO = coarseCFOTry;
        fprintf('Attempt %d: pktStart=%d coarseCFO=%.1f Hz\n', attempt, pktStart, coarseCFO);
        break;
    end

    if isempty(pktStart)
        error('Failed to lock onto valid packet.');
    end

    yPkt = y(pktStart:end);
    yPkt = frequencyOffset(yPkt, fs, -coarseCFO);

    if numel(yPkt) < pktLen
        error('Captured packet too short.');
    end

    fineCFO = wlanFineCFOEstimate(yPkt, cfgNonHT.ChannelBandwidth);
    yPkt = frequencyOffset(yPkt, fs, -fineCFO);

    symOffset = wlanSymbolTimingEstimate(yPkt(1:min(20000,end)), cfgNonHT.ChannelBandwidth);

    %% ------------ Channel estimation ------------
    rxLLTF = yPkt(ind.LLTF(1):ind.LLTF(2), :);
    lltfDemod = wlanLLTFDemodulate(rxLLTF, cfgNonHT);
    chanEst = wlanLLTFChannelEstimate(lltfDemod, cfgNonHT);
    noiseEst = wlanLLTFNoiseEstimate(lltfDemod);

    ofdmInfo = wlanNonHTOFDMInfo('NonHT-Data', cfgNonHT);
    Nfft = ofdmInfo.FFTLength;
    Ncp  = ofdmInfo.CPLength;
    activeFFT  = ofdmInfo.ActiveFFTIndices;
    pilotIdx52 = ofdmInfo.PilotIndices;
    dataIdx52  = ofdmInfo.DataIndices;

    symLen = Nfft + Ncp;

    rxData = yPkt(ind.NonHTData(1):ind.NonHTData(2), :);
    numSym = floor(size(rxData,1)/symLen);

    H = squeeze(chanEst(:,1,1));

    rxPilotsEq = zeros(numel(pilotIdx52), numSym);
    rxDataEq   = zeros(numel(dataIdx52),  numSym);

    for ell = 1:numSym
        seg = rxData((ell-1)*symLen + (1:symLen), 1);
        seg = seg(Ncp+1:end);
        Y = fft(seg, Nfft);
        Yact = Y(activeFFT);
        Zeq = Yact ./ (H + 1e-12);
        rxPilotsEq(:,ell) = Zeq(pilotIdx52);
        rxDataEq(:,ell)   = Zeq(dataIdx52);
    end

    %% ------------ Differential DPPA detector ------------
    u_meas = zeros(1,numSym);
    phiH1_deg = zeros(1,numSym);
    phiH0_deg = zeros(1,numSym);

    for ell = 1:numSym
        z = rxPilotsEq(:,ell);
        p_nom = txPilotsBase(:,ell);

        u_meas(ell) = (p_nom' * z) / (norm(p_nom)^2 + 1e-12);

        nonce_bits_H1 = nonce_with_symbol(nonce_frame, ell, 16);
        C_eff_H1      = mix_challenge_nonce(C_seed, nonce_bits_H1, length(C_seed));
        R_H1          = get_puf_bits(node_claimed, C_eff_H1, B);
        seed_H1       = hash_bits_u32(R_H1);
        phiH1_deg(ell)= phase_from_seed(seed_H1, phi_choices_deg);

        nonce_frame_H0 = [1 - nonce_frame(1), nonce_frame(2:end)];
        nonce_bits_H0  = nonce_with_symbol(nonce_frame_H0, ell, 16);
        C_eff_H0       = mix_challenge_nonce(C_seed, nonce_bits_H0, length(C_seed));
        R_H0           = get_puf_bits(node_claimed, C_eff_H0, B);
        seed_H0        = hash_bits_u32(R_H0);
        phiH0_deg(ell) = phase_from_seed(seed_H0, phi_choices_deg);
    end

    Td_H1 = zeros(1,numSym-1);
    Td_H0 = zeros(1,numSym-1);
    dErr_H1_deg = zeros(1,numSym-1);
    dErr_H0_deg = zeros(1,numSym-1);

    for ell = 2:numSym
        d_meas = u_meas(ell) * conj(u_meas(ell-1));
        dphi_meas = angle(d_meas);

        dphi_H1 = deg2rad(phiH1_deg(ell) - phiH1_deg(ell-1));
        dphi_H0 = deg2rad(phiH0_deg(ell) - phiH0_deg(ell-1));

        e1 = angle(exp(1j*(dphi_meas - dphi_H1)));
        e0 = angle(exp(1j*(dphi_meas - dphi_H0)));

        Td_H1(ell-1) = cos(e1);
        Td_H0(ell-1) = cos(e0);

        dErr_H1_deg(ell-1) = rad2deg(e1);
        dErr_H0_deg(ell-1) = rad2deg(e0);
    end

    %% ------------ DPPA frame decision ------------
    Ldiff = numel(Td_H1);
    VH1 = sum(Td_H1 >= tau_d);
    VH0 = sum(Td_H0 >= tau_d);

    acceptH1 = (VH1 >= ceil(alpha * Ldiff));
    acceptH0 = (VH0 >= ceil(alpha * Ldiff));

    %% ------------ BER ------------
    [rxPSDU, ~, ~] = wlanNonHTDataRecover(rxData, chanEst, noiseEst, cfgNonHT);
    L = min(numel(psdu), numel(rxPSDU));
    numErr = sum(rxPSDU(1:L) ~= psdu(1:L));
    ber = numErr / L;

    %% ------------ Pack results ------------
    res = struct();
    res.fc = fc;
    res.ber = ber;
    res.numErr = numErr;
    res.Lbits = L;
    res.TdH1_mean = mean(Td_H1);
    res.TdH0_mean = mean(Td_H0);
    res.TpH1_mean = mean(Td_H1);
    res.TpH0_mean = mean(Td_H0);
    res.dErrH1_mean_deg = mean(abs(dErr_H1_deg));
    res.dErrH0_mean_deg = mean(abs(dErr_H0_deg));
    res.VH1 = VH1;
    res.VH0 = VH0;
    res.Ldiff = Ldiff;
    res.acceptH1 = acceptH1;
    res.acceptH0 = acceptH0;
    res.pktStart = pktStart;
    res.coarseCFO = coarseCFO;
    res.fineCFO = fineCFO;
    res.rxMaxAmp = rxMaxAmp;
    res.Td_H1 = Td_H1;
    res.Td_H0 = Td_H0;
    res.TdH1_mean = res.TpH1_mean;
    res.TdH0_mean = res.TpH0_mean;
   
end

function S = pilotSignSequence(numSym, seed)
% Returns a 4 x numSym matrix of pilot sign patterns (+1/-1)
% Deterministic from the seed.

    rng(double(seed), 'twister');

    S = sign(randn(4, numSym));
    S(S==0) = 1;
end

function s = pilot_signs_from_seed(seed_u32, Npilots)
% pilot_signs_from_seed
% Deterministically generate Npilots signs in {+1,-1} from a uint32 seed.

    rng(double(seed_u32), 'twister');

    s = sign(randn(Npilots,1));
    s(s == 0) = 1;
end

function phi_deg = phase_from_seed(seed_u32, phi_choices_deg)
    idx_phi = 1 + mod(double(seed_u32), numel(phi_choices_deg));
    phi_deg = phi_choices_deg(idx_phi);
end

function cleanup_radios(tx, rx)
    try
        release(tx);
    catch
    end
    try
        release(rx);
    catch
    end
end