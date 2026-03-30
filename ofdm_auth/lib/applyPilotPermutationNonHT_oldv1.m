function y = applyPilotPermutationNonHT(waveform, cfgNonHT, node, C_seed, nonce_frame, B)

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

    for ell = 1:numSym
        idx = (ell-1)*symLen + (1:symLen);
        xcp = xData(idx);
        x   = xcp(Ncp+1:end);

        X64 = fft(x, Nfft);
        X52 = X64(activeFFT);

        % Current nominal pilot values for this symbol
        p_nom = X52(pilotIdx52);

        % Derive per-symbol seed from PUF pipeline
        nonce_bits = nonce_with_symbol(nonce_frame, ell, 16);
        C_eff      = mix_challenge_nonce(C_seed, nonce_bits, length(C_seed));
        R_bits     = get_puf_bits(node, C_eff, B);
        seed_u32   = hash_bits_u32(R_bits);

        % Deterministic permutation of the 4 pilot slots
        perm_idx = perm_from_seed(seed_u32, numel(pilotIdx52));

        if ell <= 3
            fprintf('TX symbol %d seed=%u perm=[%d %d %d %d]\n', ...
                ell, seed_u32, perm_idx(1), perm_idx(2), perm_idx(3), perm_idx(4));
        end

        % Permute the 4 pilot values
        X52(pilotIdx52) = p_nom(perm_idx);
        % Restricted permutation codebook: only mild pilot swaps
        % perm_codebook = [
        %     1 2 3 4;   % identity
        %     2 1 3 4;   % swap pair (1,2)
        %     1 2 4 3;   % swap pair (3,4)
        %     4 2 3 1    % swap pair (1,4) optional
        % ];
        % 
        % code_idx = 1 + mod(double(seed_u32), size(perm_codebook,1));
        % perm_idx = perm_codebook(code_idx,:);
        % 
        % if ell <= 3
        %     fprintf('TX symbol %d seed=%u perm=[%d %d %d %d]\n', ...
        %         ell, seed_u32, perm_idx(1), perm_idx(2), perm_idx(3), perm_idx(4));
        % end
        % 
        % X52(pilotIdx52) = p_nom(perm_idx);

        X64(activeFFT) = X52;
        x_mod = ifft(X64, Nfft);
        xcp_mod = [x_mod(end-Ncp+1:end); x_mod];

        xData(idx) = xcp_mod;
    end

    y(dataStart:dataEnd,1) = xData;
end