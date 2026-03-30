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
        waveform = applyPilotPermutationNonHT(waveform, cfgNonHT, node, C_seed, nonce_frame, B);
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

    %% ------------ Permutation detector ------------
    Tp_H1 = zeros(1,numSym);
    Tp_H0 = zeros(1,numSym);

    perm_codebook = [
    1 2 3 4;
    2 1 3 4;
    1 2 4 3;
    4 2 3 1
    ];

    for ell = 1:numSym
        z = rxPilotsEq(:,ell);
        p_nom = txPilotsBase(:,ell);

        % H1 predicted permutation
        nonce_bits_H1 = nonce_with_symbol(nonce_frame, ell, 16);
        C_eff_H1      = mix_challenge_nonce(C_seed, nonce_bits_H1, length(C_seed));
        R_H1          = get_puf_bits(node_claimed, C_eff_H1, B);
        seed_H1       = hash_bits_u32(R_H1);
        perm_H1       = perm_from_seed(seed_H1, numel(pilotIdx52));
        %idx_H1 = 1 + mod(double(seed_H1), size(perm_codebook,1));
        %perm_H1 = perm_codebook(idx_H1,:);
        p_H1          = p_nom(perm_H1);

        % H0 wrong/stale nonce hypothesis
        nonce_frame_H0 = [1 - nonce_frame(1), nonce_frame(2:end)];
        nonce_bits_H0  = nonce_with_symbol(nonce_frame_H0, ell, 16);
        C_eff_H0       = mix_challenge_nonce(C_seed, nonce_bits_H0, length(C_seed));
        R_H0           = get_puf_bits(node_claimed, C_eff_H0, B);
        seed_H0        = hash_bits_u32(R_H0);
        perm_H0        = perm_from_seed(seed_H0, numel(pilotIdx52));
        %idx_H0 = 1 + mod(double(seed_H0), size(perm_codebook,1));
        %perm_H0 = perm_codebook(idx_H0,:);
        p_H0           = p_nom(perm_H0);

        Tp_H1(ell) = real((p_H1' * z) / (norm(p_H1)^2 + 1e-12));
        Tp_H0(ell) = real((p_H0' * z) / (norm(p_H0)^2 + 1e-12));
    end

    %% ------------ Frame decision ------------
    VH1 = sum(Tp_H1 >= tau_p);
    VH0 = sum(Tp_H0 >= tau_p);
    acceptH1 = (VH1 >= ceil(alpha * numSym));
    acceptH0 = (VH0 >= ceil(alpha * numSym));

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
    res.TpH1_mean = mean(Tp_H1);
    res.TpH0_mean = mean(Tp_H0);
    res.VH1 = VH1;
    res.VH0 = VH0;
    res.Lsym = numSym;
    res.acceptH1 = acceptH1;
    res.acceptH0 = acceptH0;
    res.pktStart = pktStart;
    res.coarseCFO = coarseCFO;
    res.fineCFO = fineCFO;
    res.rxMaxAmp = rxMaxAmp;
    res.Tp_H1 = Tp_H1;
    res.Tp_H0 = Tp_H0;
end