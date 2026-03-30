function y = applyPilotPermutationNonHT(waveform, cfgNonHT, node, C_seed, nonce_frame, B, perm_cfg)
% applyPilotPermutationNonHT
% Constrained payload-pilot permutation for OTA WLAN.
%
% perm_cfg fields:
%   fixedPilotLocalIdx : pilot local indices left unchanged
%   authPilotLocalIdx  : pilot local indices used for auth permutation

    y = waveform;

    ind = wlanFieldIndices(cfgNonHT);

    ofdmInfo   = wlanNonHTOFDMInfo('NonHT-Data', cfgNonHT);
    Nfft       = ofdmInfo.FFTLength;
    Ncp        = ofdmInfo.CPLength;
    symLen     = Nfft + Ncp;
    activeFFT  = ofdmInfo.ActiveFFTIndices;
    pilotIdx52 = ofdmInfo.PilotIndices;

    dataStart = ind.NonHTData(1);
    dataEnd   = ind.NonHTData(2);

    xData = y(dataStart:dataEnd, 1);

    if mod(numel(xData), symLen) ~= 0
        error('NonHT-Data length is not an integer multiple of OFDM symbol length.');
    end

    numSym = numel(xData) / symLen;

    fixedLocal = perm_cfg.fixedPilotLocalIdx(:).';
    authLocal  = perm_cfg.authPilotLocalIdx(:).';

    assert(all(ismember([fixedLocal authLocal], 1:numel(pilotIdx52))), ...
        'Pilot local indices must be within the 4 WLAN pilot positions.');
    assert(isempty(intersect(fixedLocal, authLocal)), ...
        'fixedPilotLocalIdx and authPilotLocalIdx must be disjoint.');
    %assert(numel(authLocal) == 2, 'This constrained design assumes exactly 2 auth pilots.');
    assert(numel(authLocal) >= 2, 'Need at least 2 auth pilots.');

    for ell = 1:numSym
        idx = (ell-1)*symLen + (1:symLen);
        xcp = xData(idx);
        x   = xcp(Ncp+1:end);

        X64 = fft(x, Nfft);
        X52 = X64(activeFFT);

        p_nom = X52(pilotIdx52);   % 4x1 nominal pilot vector
        p_mod = p_nom;             % start from nominal

        % ----------- Inject pilot asymmetry (CRITICAL) -----------
        % These weights create distinguishable pilot identities
        pilotWeights = [1.0, 0.35, 0.8, 0.2];   % tuneable
        
        p_nom = p_nom .* pilotWeights(:);   % apply asymmetry
        p_mod = p_nom;

        % Derive per-symbol seed
        nonce_bits = nonce_with_symbol(nonce_frame, ell, 16);
        C_eff      = mix_challenge_nonce(C_seed, nonce_bits, length(C_seed));
        R_bits     = get_puf_bits(node, C_eff, B);
        seed_u32   = hash_bits_u32(R_bits);

        % 2-state swap only on the auth subset
        %perm_auth = perm_state_from_seed(seed_u32);
        perm_auth = perm_from_seed(seed_u32, numel(authLocal));


        if ell <= 3
            fprintf('TX symbol %d seed=%u permAuth=[%d %d]\n', ...
                ell, seed_u32, perm_auth(1), perm_auth(2));
        end

        %p_mod(authLocal) = p_nom(authLocal(perm_auth));
        p_mod(authLocal) = p_nom(authLocal(perm_auth));

        X52(pilotIdx52) = p_mod;
        X64(activeFFT) = X52;

        x_mod = ifft(X64, Nfft);
        xcp_mod = [x_mod(end-Ncp+1:end); x_mod];
        xData(idx) = xcp_mod;
    end

    y(dataStart:dataEnd,1) = xData;
end