function txPilotsBase = extractPayloadPilotsFromWaveform(waveform, cfgNonHT)
% extractPayloadPilotsFromWaveform
% Extract nominal payload pilot values from an uncoded Non-HT WLAN waveform.
%
% Inputs:
%   waveform  : output of wlanWaveformGenerator(psdu, cfgNonHT, ...)
%   cfgNonHT  : wlanNonHTConfig object
%
% Output:
%   txPilotsBase : [Npilots x Nsym] nominal pilot matrix from the payload

    ind = wlanFieldIndices(cfgNonHT);
    ofdmInfo = wlanNonHTOFDMInfo('NonHT-Data', cfgNonHT);

    Nfft       = ofdmInfo.FFTLength;
    Ncp        = ofdmInfo.CPLength;
    symLen     = Nfft + Ncp;
    activeFFT  = ofdmInfo.ActiveFFTIndices;
    pilotIdx52 = ofdmInfo.PilotIndices;

    % Extract Non-HT data field only
    xData = waveform(ind.NonHTData(1):ind.NonHTData(2), 1);

    if mod(numel(xData), symLen) ~= 0
        error('NonHT-Data length is not an integer multiple of OFDM symbol length.');
    end

    numSym = numel(xData) / symLen;
    txPilotsBase = zeros(numel(pilotIdx52), numSym);

    for ell = 1:numSym
        idx = (ell-1)*symLen + (1:symLen);
        xcp = xData(idx);
        x   = xcp(Ncp+1:end);   % remove CP

        X64 = fft(x, Nfft);
        X52 = X64(activeFFT);

        txPilotsBase(:, ell) = X52(pilotIdx52);
    end
end