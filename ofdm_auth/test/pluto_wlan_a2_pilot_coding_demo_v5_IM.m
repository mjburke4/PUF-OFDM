clear; clc; %close all;
addpath 'C:\Users\Mike B\OneDrive\Desktop\Grad School\Phd\Experiments\PUF_OFDM\PUF OFDM Repository\ofdm_auth\lib'

Ntrials = 10;
fc_list = [350e6];
tau_im  = 0.20;     % starter threshold, tune later
alpha   = 0.75;


fprintf('\n===== STARTING IM OTA RUN =====\n');

for k = 1:Ntrials
    fc_this = fc_list(mod(k-1, numel(fc_list)) + 1);
    fprintf('\nTrial %d / %d  fc = %.3f MHz\n', k, Ntrials, fc_this/1e6);

    try
        res = run_one_im_trial(fc_this, tau_im, alpha);
        fprintf('BER=%.3e, Tim_H1=%.3f, Tim_H0=%.3f, acc=%d/%d\n', ...
            res.ber, res.TimH1_mean, res.TimH0_mean, res.acceptH1, res.acceptH0);
    catch ME
        fprintf(2, 'Trial %d FAILED: %s\n', k, ME.message);
    end
end


Ntrials = 20;
fc_list = [350e6];
tau_im  = 0.20;      % starter threshold, tune later
alpha   = 0.75;

results = repmat(struct( ...
    'trial', [], ...
    'fc', [], ...
    'ber', [], ...
    'numErr', [], ...
    'Lbits', [], ...
    'TimH1_mean', [], ...
    'TimH0_mean', [], ...
    'VH1', [], ...
    'VH0', [], ...
    'Lsym', [], ...
    'Lauth', [], ...   % <-- ADD THIS
    'acceptH1', [], ...
    'acceptH0', [], ...
    'pktStart', [], ...
    'coarseCFO', [], ...
    'fineCFO', [], ...
    'rxMaxAmp', [], ...
    'Tim_H1', [], ...
    'Tim_H0', [], ...
    'ok', false, ...
    'errmsg', ''), Ntrials, 1);

fprintf('\n===== STARTING IM OTA RUN =====\n');

for k = 1:Ntrials
    fc_this = fc_list(mod(k-1, numel(fc_list)) + 1);

    fprintf('\n----------------------------------------\n');
    fprintf('Trial %d / %d   fc = %.3f MHz\n', k, Ntrials, fc_this/1e6);
    fprintf('----------------------------------------\n');

    try
        res = run_one_im_trial(fc_this, tau_im, alpha);

        results(k).trial      = k;
        results(k).fc         = res.fc;
        results(k).ber        = res.ber;
        results(k).numErr     = res.numErr;
        results(k).Lbits      = res.Lbits;
        results(k).TimH1_mean = res.TimH1_mean;
        results(k).TimH0_mean = res.TimH0_mean;
        results(k).VH1        = res.VH1;
        results(k).VH0        = res.VH0;
        results(k).Lsym       = res.Lsym;
        results(k).Lauth = res.Lauth;
        results(k).acceptH1   = res.acceptH1;
        results(k).acceptH0   = res.acceptH0;
        results(k).pktStart   = res.pktStart;
        results(k).coarseCFO  = res.coarseCFO;
        results(k).fineCFO    = res.fineCFO;
        results(k).rxMaxAmp   = res.rxMaxAmp;
        results(k).Tim_H1     = res.Tim_H1;
        results(k).Tim_H0     = res.Tim_H0;
        results(k).ok         = true;
        results(k).errmsg     = '';

        fprintf('Trial %d summary: BER=%.3e, Tim_H1=%.3f, Tim_H0=%.3f, acc(H1/H0)=%d/%d\n', ...
            k, res.ber, res.TimH1_mean, res.TimH0_mean, res.acceptH1, res.acceptH0);

    catch ME
        results(k).trial = k;
        results(k).fc    = fc_this;
        results(k).ok    = false;
        results(k).errmsg = ME.message;

        fprintf(2, 'Trial %d FAILED: %s\n', k, ME.message);
    end
end

fprintf('\n===== MULTI-TRIAL RUN COMPLETE =====\n');

ok_idx = find([results.ok]);

if isempty(ok_idx)
    error('No successful trials completed.');
end

ber_vec   = [results(ok_idx).ber];
TimH1_vec = [results(ok_idx).TimH1_mean];
TimH0_vec = [results(ok_idx).TimH0_mean];
accH1_vec = [results(ok_idx).acceptH1];
accH0_vec = [results(ok_idx).acceptH0];
Lsym_vec  = [results(ok_idx).Lsym];
Lauth_vec = [results(ok_idx).Lauth];

fprintf('\n===== SUMMARY OVER %d SUCCESSFUL TRIALS =====\n', numel(ok_idx));
fprintf('Mean BER              = %.3e\n', mean(ber_vec));
fprintf('Median BER            = %.3e\n', median(ber_vec));
fprintf('Mean Tim_H1           = %.3f\n', mean(TimH1_vec));
fprintf('Mean Tim_H0           = %.3f\n', mean(TimH0_vec));
fprintf('H1 acceptance rate    = %.3f\n', mean(accH1_vec));
fprintf('H0 acceptance rate    = %.3f\n', mean(accH0_vec));
fprintf('Mean payload Lsym     = %.1f\n', mean(Lsym_vec));
fprintf('Mean auth-window Lauth = %.1f\n', mean(Lauth_vec));

T = table( ...
    [results.trial].', ...
    [results.fc].'/1e6, ...
    [results.ok].', ...
    [results.ber].', ...
    [results.TimH1_mean].', ...
    [results.TimH0_mean].', ...
    [results.acceptH1].', ...
    [results.acceptH0].', ...
    [results.Lsym].', ...
    'VariableNames', {'trial','fc_MHz','ok','BER','TimH1_mean','TimH0_mean','acceptH1','acceptH0','Lsym'});

disp(T);

figure;
plot(ok_idx, ber_vec, 'o-','LineWidth',1.5);
grid on;
xlabel('Trial index');
ylabel('BER');
title('BER across IM OTA trials');

figure;
plot(ok_idx, TimH1_vec, 'o-','LineWidth',1.5); hold on;
plot(ok_idx, TimH0_vec, 's-','LineWidth',1.5);
grid on;
legend('Mean T_{IM,H1}','Mean T_{IM,H0}','Location','best');
xlabel('Trial index');
ylabel('Mean IM score');
title('IM score separation across OTA trials');

function res = run_one_im_trial(fc, tau_im, alpha)

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

    %% ------------ PUF setup ------------
    U = 1;
    u0 = 1;
    m = 64;
    B = 64;

    nodes = init_node_population(U, m, 'puf_bits', 128, 'seed', 7);
    node = nodes(u0);
    node_claimed = node;

    C_seed = randi([0 1], 1, node.m);
    nonce_frame = randi([0 1], 1, 32);

    %% ------------ IM configuration ------------
    im_cfg.fixedPilotLocalIdx = [1 4];   % leave these unchanged for WLAN tracking
    im_cfg.imPilotLocalIdx    = [2 3];   % only these carry IM auth
    im_cfg.K_active_im        = 1;       % one of the two is "active"
    im_cfg.inactiveScale      = 0.25;    % gentle suppression first; later can try 0
    %im_cfg.inactiveScale = 0.5;

    %% ------------ Waveform generation ------------
    waveform_uncoded = wlanWaveformGenerator(psdu, cfgNonHT, 'IdleTime', idleTime);
    txPilotsBase = extractPayloadPilotsFromWaveform(waveform_uncoded, cfgNonHT);

    waveform = waveform_uncoded;
    if enablePilotCoding
        waveform = applyPilotIMNonHT(waveform, cfgNonHT, node, C_seed, nonce_frame, B, im_cfg);
    end

    guard = zeros(round(500e-6*fs),1);
    waveform = [guard; waveform; guard];
    waveform = 0.6 * waveform / max(abs(waveform)+1e-12);

    %% ------------ Pluto TX/RX ------------
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

    %% ------------ TX/RX capture with retry ------------
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
        break
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

    %% ------------ Channel estimate ------------
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

    %% ------------ IM detector ------------
    Tim_H1 = zeros(1,numSym);
    Tim_H0 = zeros(1,numSym);

    fixedLocal = im_cfg.fixedPilotLocalIdx(:).';
    imLocal    = im_cfg.imPilotLocalIdx(:).';
    K_active   = im_cfg.K_active_im;

    for ell = 1:numSym
        z = rxPilotsEq(:,ell);
        e = abs(z).^2;

        % H1 predicted active IM subset
        nonce_bits_H1 = nonce_with_symbol(nonce_frame, ell, 16);
        C_eff_H1      = mix_challenge_nonce(C_seed, nonce_bits_H1, length(C_seed));
        R_H1          = get_puf_bits(node_claimed, C_eff_H1, B);
        seed_H1       = hash_bits_u32(R_H1);
        mask_H1       = mask_from_seed_subset(seed_H1, numel(imLocal), K_active);

        activeH1 = imLocal(mask_H1 > 0);
        inactiveH1 = imLocal(mask_H1 == 0);

        Ea1 = mean(e(activeH1));
        Ei1 = mean(e(inactiveH1));
        Tim_H1(ell) = (Ea1 - Ei1) / (Ea1 + Ei1 + 1e-12);

        % H0 wrong/stale nonce hypothesis
        nonce_frame_H0 = [1 - nonce_frame(1), nonce_frame(2:end)];
        nonce_bits_H0  = nonce_with_symbol(nonce_frame_H0, ell, 16);
        C_eff_H0       = mix_challenge_nonce(C_seed, nonce_bits_H0, length(C_seed));
        R_H0           = get_puf_bits(node_claimed, C_eff_H0, B);
        seed_H0        = hash_bits_u32(R_H0);
        mask_H0        = mask_from_seed_subset(seed_H0, numel(imLocal), K_active);

        activeH0 = imLocal(mask_H0 > 0);
        inactiveH0 = imLocal(mask_H0 == 0);

        Ea0 = mean(e(activeH0));
        Ei0 = mean(e(inactiveH0));
        Tim_H0(ell) = (Ea0 - Ei0) / (Ea0 + Ei0 + 1e-12);
    end

    % ------------ Authentication window length ------------
    Lauth = 16;   % use first 16 payload OFDM symbols for auth decision
    Lauth_eff = min(Lauth, numSym);
    
    Tim_H1_auth = Tim_H1(1:Lauth_eff);
    Tim_H0_auth = Tim_H0(1:Lauth_eff);
    
    VH1 = sum(Tim_H1_auth >= tau_im);
    VH0 = sum(Tim_H0_auth >= tau_im);
    
    acceptH1 = (VH1 >= ceil(alpha * Lauth_eff));
    acceptH0 = (VH0 >= ceil(alpha * Lauth_eff));

    res.TimH1_mean = mean(Tim_H1_auth);
    res.TimH0_mean = mean(Tim_H0_auth);
    
    res.Lsym  = numSym;      % full payload symbols
    res.Lauth = Lauth_eff;   % auth symbols actually used

    %% ------------ Frame decision ------------
    % VH1 = sum(Tim_H1 >= tau_im);
    % VH0 = sum(Tim_H0 >= tau_im);
    % acceptH1 = (VH1 >= ceil(alpha * numSym));
    % acceptH0 = (VH0 >= ceil(alpha * numSym));

    %% ------------ BER ------------
    [rxPSDU, ~, ~] = wlanNonHTDataRecover(rxData, chanEst, noiseEst, cfgNonHT);
    L = min(numel(psdu), numel(rxPSDU));
    numErr = sum(rxPSDU(1:L) ~= psdu(1:L));
    ber = numErr / L;

    %% ------------ Results ------------
    res = struct();
    res.fc = fc;
    res.ber = ber;
    res.numErr = numErr;
    res.Lbits = L;
    % res.TimH1_mean = mean(Tim_H1);
    % res.TimH0_mean = mean(Tim_H0);
    res.TimH1_mean = mean(Tim_H1_auth);
    res.TimH0_mean = mean(Tim_H0_auth);
    res.VH1 = VH1;
    res.VH0 = VH0;
    res.Lsym = numSym;
    res.acceptH1 = acceptH1;
    res.acceptH0 = acceptH0;
    res.pktStart = pktStart;
    res.coarseCFO = coarseCFO;
    res.fineCFO = fineCFO;
    res.rxMaxAmp = rxMaxAmp;
    res.Tim_H1 = Tim_H1;
    res.Tim_H0 = Tim_H0;
    res.Lauth = Lauth_eff;
    res.Tim_H1_auth = Tim_H1_auth;
    res.Tim_H0_auth = Tim_H0_auth;
end