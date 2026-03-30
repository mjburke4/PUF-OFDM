function res = run_one_perm_trial(fc, tau_p, alpha)

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
    cfgNonHT.MCS = 1;  % use 2 for QPSK - Le'te do BPSK from now on though
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

    %% ------------ Constrained permutation config ------------
    perm_cfg.fixedPilotLocalIdx = [1 4];
    perm_cfg.authPilotLocalIdx  = [2 3];

    %perm_cfg.fixedPilotLocalIdx = [];        % keep only 1 fixed
    %perm_cfg.authPilotLocalIdx  = [1 2 3 4];    % 3 auth pilots

    Lauth = 64; %16;   % auth window for fairer comparison to sim

    %% ------------ Waveform generation ------------
    waveform_uncoded = wlanWaveformGenerator(psdu, cfgNonHT, 'IdleTime', idleTime);
    txPilotsBase = extractPayloadPilotsFromWaveform(waveform_uncoded, cfgNonHT);

    waveform = waveform_uncoded;
    if enablePilotCoding
        waveform = applyPilotPermutationNonHT(waveform, cfgNonHT, node, C_seed, nonce_frame, B, perm_cfg);
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

    %% ------------ TX/RX capture with retry ------------
    transmitRepeat(tx, waveform);
    pause(pauseTX);

    maxTries = 8;
    pktStart = [];
    coarseCFO = NaN;
    rxMaxAmp = NaN;
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

    %% ------------ Constrained permutation detector ------------
    Tp_H1 = zeros(1,numSym);
    Tp_H0 = zeros(1,numSym);

    fixedLocal = perm_cfg.fixedPilotLocalIdx(:).';
    authLocal  = perm_cfg.authPilotLocalIdx(:).';

    pilotWeights = [1.0, 0.35, 0.8, 0.2];   % MUST match TX

for ell = 1:numSym
    % Received equalized pilots for this symbol
    z = rxPilotsEq(:,ell);

    % Nominal TX pilot template for this symbol, with same asymmetry as TX
    p_nom = txPilotsBase(:,ell) .* pilotWeights(:);

    % -------------------------------------------------
    % Step 1: estimate common phase error from fixed pilots
    % -------------------------------------------------
    z_fix = z(fixedLocal);
    p_fix = p_nom(fixedLocal);

    theta_hat = angle(p_fix' * z_fix);

    % De-rotate all received pilots using fixed-pilot phase estimate
    z_corr = z * exp(-1j * theta_hat);

    % Auth subset after phase correction
    z_auth = z_corr(authLocal);
    p_nom_auth = p_nom(authLocal);

    % -------------------------------------------------
    % Step 2: H1 predicted auth permutation
    % -------------------------------------------------
    nonce_bits_H1 = nonce_with_symbol(nonce_frame, ell, 16);
    C_eff_H1      = mix_challenge_nonce(C_seed, nonce_bits_H1, length(C_seed));
    R_H1          = get_puf_bits(node_claimed, C_eff_H1, B);
    seed_H1       = hash_bits_u32(R_H1);

    perm_auth_H1 = perm_from_seed(seed_H1, numel(authLocal));
    p_H1_auth    = p_nom_auth(perm_auth_H1);

    % -------------------------------------------------
    % Step 3: H0 predicted auth permutation (stale/wrong nonce)
    % -------------------------------------------------
    nonce_frame_H0 = [1 - nonce_frame(1), nonce_frame(2:end)];
    nonce_bits_H0  = nonce_with_symbol(nonce_frame_H0, ell, 16);
    C_eff_H0       = mix_challenge_nonce(C_seed, nonce_bits_H0, length(C_seed));
    R_H0           = get_puf_bits(node_claimed, C_eff_H0, B);
    seed_H0        = hash_bits_u32(R_H0);

    perm_auth_H0 = perm_from_seed(seed_H0, numel(authLocal));
    p_H0_auth    = p_nom_auth(perm_auth_H0);

    % -------------------------------------------------
    % Step 4: correlator scores
    % -------------------------------------------------
    % Tp_H1(ell) = real((p_H1_auth' * z_auth) / (norm(p_H1_auth)^2 + 1e-12));
    % Tp_H0(ell) = real((p_H0_auth' * z_auth) / (norm(p_H0_auth)^2 + 1e-12));

    z_auth_n = z_auth / (norm(z_auth) + 1e-12);
    p_H1_n   = p_H1_auth / (norm(p_H1_auth) + 1e-12);
    p_H0_n   = p_H0_auth / (norm(p_H0_auth) + 1e-12);
    
    Tp_H1(ell) = real(p_H1_n' * z_auth_n);
    Tp_H0(ell) = real(p_H0_n' * z_auth_n);
    
    if ell <= 3
    fprintf('ell=%d theta_hat=%.2f deg  Tp_H1=%.3f  Tp_H0=%.3f\n', ...
        ell, rad2deg(theta_hat), Tp_H1(ell), Tp_H0(ell));
    end
end
    % for ell = 1:numSym
    %     z = rxPilotsEq(:,ell);
    %     % z_auth = z(authLocal);
    % 
    %     %p_nom = txPilotsBase(:,ell);
    %     % Test code to create weights for the pilots
    %     pilotWeights = [1.0, 0.35, 0.8, 0.2];   % MUST match TX
    %     p_nom = txPilotsBase(:,ell) .* pilotWeights(:);
    %     % 
    %     % p_nom_auth = p_nom(authLocal);
    % 
    %     % --- Use fixed pilots to estimate common phase error ---
    % 
    %     % Fixed-pilot common phase correction
    %     z_fix = z(fixedLocal);
    %     p_fix = p_nom(fixedLocal);
    % 
    %     theta_hat = angle(p_fix' * z_fix);
    %     z_corr = z * exp(-1j * theta_hat);
    % 
    %     z_auth = z_corr(authLocal);
    %     p_nom_auth = p_nom(authLocal);
    % 
    %     %theta_hat = angle(p_fix' * z_fix);   % common phase estimate from fixed pilots
    % 
    %     % De-rotate all received pilots
    %     %z_corr = z * exp(-1j * theta_hat);
    % 
    %     % Auth subset after phase correction
    %     % z_auth = z_corr(authLocal);
    %     % p_nom_auth = p_nom(authLocal);
    % 
    %     % H1 predicted auth-swap state
    %     nonce_bits_H1 = nonce_with_symbol(nonce_frame, ell, 16);
    %     C_eff_H1      = mix_challenge_nonce(C_seed, nonce_bits_H1, length(C_seed));
    %     R_H1          = get_puf_bits(node_claimed, C_eff_H1, B);
    %     seed_H1       = hash_bits_u32(R_H1);
    %     %perm_auth_H1  = perm_state_from_seed(seed_H1);
    %     %p_H1_auth     = p_nom_auth(perm_auth_H1);
    %     perm_auth_H1 = perm_from_seed(seed_H1, numel(authLocal));
    %     p_H1_auth = p_nom_auth(perm_auth_H1);
    % 
    %     % H0 wrong/stale nonce hypothesis
    %     nonce_frame_H0 = [1 - nonce_frame(1), nonce_frame(2:end)];
    %     nonce_bits_H0  = nonce_with_symbol(nonce_frame_H0, ell, 16);
    %     C_eff_H0       = mix_challenge_nonce(C_seed, nonce_bits_H0, length(C_seed));
    %     R_H0           = get_puf_bits(node_claimed, C_eff_H0, B);
    %     seed_H0        = hash_bits_u32(R_H0);
    %     %perm_auth_H0   = perm_state_from_seed(seed_H0);
    %     %p_H0_auth      = p_nom_auth(perm_auth_H0);
    %     perm_auth_H0 = perm_from_seed(seed_H0, numel(authLocal));
    %     p_H0_auth = p_nom_auth(perm_auth_H0);
    % 
    %     Tp_H1(ell) = real((p_H1_auth' * z_auth) / (norm(p_H1_auth)^2 + 1e-12));
    %     Tp_H0(ell) = real((p_H0_auth' * z_auth) / (norm(p_H0_auth)^2 + 1e-12));
    % end

    %% ------------ Authentication window ------------
    Lauth_eff = min(Lauth, numSym);

    Tp_H1_auth = Tp_H1(1:Lauth_eff);
    Tp_H0_auth = Tp_H0(1:Lauth_eff);

    VH1 = sum(Tp_H1_auth >= tau_p);
    VH0 = sum(Tp_H0_auth >= tau_p);

    acceptH1 = (VH1 >= ceil(alpha * Lauth_eff));
    acceptH0 = (VH0 >= ceil(alpha * Lauth_eff));

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
    res.TpH1_mean = mean(Tp_H1_auth);
    res.TpH0_mean = mean(Tp_H0_auth);
    res.VH1 = VH1;
    res.VH0 = VH0;
    res.Lsym = numSym;
    res.Lauth = Lauth_eff;
    res.acceptH1 = acceptH1;
    res.acceptH0 = acceptH0;
    res.pktStart = pktStart;
    res.coarseCFO = coarseCFO;
    res.fineCFO = fineCFO;
    res.rxMaxAmp = rxMaxAmp;
    res.Tp_H1 = Tp_H1;
    res.Tp_H0 = Tp_H0;
    res.Tp_H1_auth = Tp_H1_auth;
    res.Tp_H0_auth = Tp_H0_auth;

    % temporary compatibility aliases if old summary code still expects Td names
    res.TdH1_mean = res.TpH1_mean;
    res.TdH0_mean = res.TpH0_mean;
end