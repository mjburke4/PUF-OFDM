function y = applyPilotIMNonHT(waveform, cfgNonHT, node, C_seed, nonce_frame, B, im_cfg)
% applyPilotIMNonHT
% Apply conservative pilot-domain IM to payload OFDM symbols.
%
% im_cfg fields:
%   fixedPilotLocalIdx   : local pilot indices always left unchanged
%   imPilotLocalIdx      : local pilot indices controlled by IM
%   K_active_im          : how many IM pilots remain active each symbol
%   inactiveScale        : scale factor for "inactive" pilots (0 for hard-off,
%                          or small nonzero like 0.2 for gentler suppression)

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

    fixedLocal = im_cfg.fixedPilotLocalIdx(:).';
    imLocal    = im_cfg.imPilotLocalIdx(:).';
    K_active   = im_cfg.K_active_im;
    inactiveScale = im_cfg.inactiveScale;

    assert(all(ismember([fixedLocal imLocal], 1:numel(pilotIdx52))), ...
        'Pilot local indices must be within the 4 WLAN pilot positions.');
    assert(isempty(intersect(fixedLocal, imLocal)), ...
        'fixedPilotLocalIdx and imPilotLocalIdx must be disjoint.');
    assert(K_active <= numel(imLocal), 'K_active_im must be <= number of IM pilots.');

    for ell = 1:numSym
        idx = (ell-1)*symLen + (1:symLen);
        xcp = xData(idx);
        x   = xcp(Ncp+1:end);

        X64 = fft(x, Nfft);
        X52 = X64(activeFFT);

        p_nom = X52(pilotIdx52);

        % Derive per-symbol seed
        nonce_bits = nonce_with_symbol(nonce_frame, ell, 16);
        C_eff      = mix_challenge_nonce(C_seed, nonce_bits, length(C_seed));
        R_bits     = get_puf_bits(node, C_eff, B);
        seed_u32   = hash_bits_u32(R_bits);

        % Build a K-of-N mask only on the IM subset
        mask_im = mask_from_seed_subset(seed_u32, numel(imLocal), K_active);

        if ell <= 3
            fprintf('TX symbol %d seed=%u IMmask=[%s]\n', ...
                ell, seed_u32, num2str(mask_im));
        end

        % Start with nominal pilots
        p_mod = p_nom;

        % Fixed pilots unchanged
        % IM pilots: active stay same, inactive attenuated or zeroed
        for ii = 1:numel(imLocal)
            loc = imLocal(ii);
            if mask_im(ii) == 0
                p_mod(loc) = inactiveScale * p_nom(loc);
            end
        end

        X52(pilotIdx52) = p_mod;
        X64(activeFFT) = X52;

        x_mod = ifft(X64, Nfft);
        xcp_mod = [x_mod(end-Ncp+1:end); x_mod];
        xData(idx) = xcp_mod;
    end

    y(dataStart:dataEnd,1) = xData;
end