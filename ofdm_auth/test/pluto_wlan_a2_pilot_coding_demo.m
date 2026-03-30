%% pluto_wlan_a2_pilot_coding_demo.m
% A2-0 baseline non-HT WLAN OTA with Pluto
% A2-1 pilot sign coding in DATA symbols (safe modification)
%
% Requires WLAN Toolbox + Pluto support package.

clear; clc; close all;

%% ------------ USER SETTINGS ------------
txID = 'usb:0';
rxID = 'usb:1';
%configurePlutoRadio('AD9364');

fc = 940e6;          % recommended first (cleaner RF). Change to 2.462e9 later.

fs = 20e6;           % non-HT 20 MHz baseband
txGain = -3;
rxGain = 45;

numPkts = 3;         % how many packets to capture and process
enablePilotCoding = true;

% "Seed" for pilot sign coding (deterministic)
seed = uint32(12345);

% Pluto capture settings
samplesPerFrame = 300000;  % enough to catch a packet burst
pauseTX = 0.2;
%% --------------------------------------

%% ------------ WLAN non-HT config ------------
cfgNonHT = wlanNonHTConfig;
cfgNonHT.ChannelBandwidth = 'CBW20';
cfgNonHT.MCS = 1;                     
cfgNonHT.NumTransmitAntennas = 1;

% payload
psduLen = 400;                         % bytes

cfgNonHT.PSDULength = psduLen;   % bytes

%psdu = randi([0 1], 8*psduLen, 1, 'int8');
psdu = randi([0 1], 8*cfgNonHT.PSDULength, 1, 'int8');

%ind = wlanFieldIndices(cfgNonHT);
ind = wlanFieldIndices(cfgNonHT);

pktLen = ind.NonHTData(2);   % total samples from packet start through end of data field
fprintf('Required packet length after pktStart: %d samples\n', pktLen);

% Generate reference waveform
idleTime = 200e-6;                      % idle between packets
waveform = wlanWaveformGenerator(psdu, cfgNonHT, ...
    'IdleTime', idleTime);

% Optional: apply pilot sign coding to DATA symbols (A2-1)
if enablePilotCoding
    waveform = applyPilotSignCodingNonHT(waveform, cfgNonHT, seed);
end

% Normalize to avoid clipping
guard = zeros(round(500e-6*fs),1);  % 500 us of silence at sample rate fs
waveform = [guard; waveform; guard];
waveform = 0.6 * waveform / max(abs(waveform)+1e-12);

fprintf('Waveform length: %d samples\n', numel(waveform));

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

%% ------------ Transmit repeat and capture ------------
transmitRepeat(tx, waveform);
pause(pauseTX);

fprintf('Capturing...\n');
% ---------------- Robust capture with retry ----------------
maxTries = 8;
for attempt = 1:maxTries

    y = rx();

    fprintf('RX max amplitude %.3f\n', max(abs(y)));

    pktStart = wlanPacketDetect(y, cfgNonHT.ChannelBandwidth);

    if isempty(pktStart)
        fprintf('Attempt %d: no packet detected\n', attempt);
        continue
    end

    pktStart = pktStart(1);

    yPkt = y(pktStart:end);

    coarseCFO = wlanCoarseCFOEstimate(yPkt, cfgNonHT.ChannelBandwidth);

    % Reject obviously bad detections
    if abs(coarseCFO) > 50e3
        fprintf('Attempt %d: bad packet detect (CFO %.1f Hz)\n', attempt, coarseCFO);
        continue
    end

    fprintf('Attempt %d: pktStart=%d coarseCFO=%.1f Hz\n', attempt, pktStart, coarseCFO);

    break
end

if attempt == maxTries
    error('Failed to lock onto valid packet after %d tries', maxTries);
end

% Apply coarse CFO correction
yPkt = frequencyOffset(yPkt, fs, -coarseCFO);

% y = rx();
% 
% release(tx);
% release(rx);
% 
% fprintf('RX max amplitude %.3f\n', max(abs(y)));
% 
% %% ------------ Receiver chain (based on MATLAB example blocks) ------------
% % Packet detect (L-STF)
% pktStart = wlanPacketDetect(y, cfgNonHT.ChannelBandwidth);
% 
% if isempty(pktStart)
%     error('No packet detected. Increase RX gain or move radios closer.');
% end
% 
% pktStart = pktStart(1);  % use first packet
% fprintf('Packet start index: %d\n', pktStart);
% 
% % Extract from pktStart onward
% 
% % Only take enough samples for one packet plus margin
% ind = wlanFieldIndices(cfgNonHT);
% pktLen = ind.NonHTData(2);
% 
% % After pktStart is found:
% margin = 100000; %60000;   % was 20000
% yPkt = y(pktStart : min(pktStart + pktLen + margin - 1, numel(y)));

if numel(yPkt) < pktLen
    error('Captured packet too short: need %d samples after pktStart, got %d. Increase SamplesPerFrame.', ...
        pktLen, numel(yPkt));
end

% Coarse CFO estimate/correct using L-STF
coarseCFO = wlanCoarseCFOEstimate(yPkt, cfgNonHT.ChannelBandwidth);
yPkt = frequencyOffset(yPkt, fs, -coarseCFO);

% Symbol timing using L-LTF
% timingSearch = yPkt(1:20000);   % only search first 20k samples
% symOffset = wlanSymbolTimingEstimate(timingSearch, cfgNonHT.ChannelBandwidth);

timingSearchLen = min(20000, numel(yPkt));
timingSearch = yPkt(1:timingSearchLen);

symOffset = wlanSymbolTimingEstimate(timingSearch, cfgNonHT.ChannelBandwidth);

% After timing correction, ensure enough samples remain
if numel(yPkt) < pktLen
    error('After symOffset correction, packet is too short. Need %d, have %d. Increase margin.', ...
        pktLen, numel(yPkt));
end

% Fine CFO estimate/correct using L-LTF
fineCFO = wlanFineCFOEstimate(yPkt, cfgNonHT.ChannelBandwidth);
yPkt = frequencyOffset(yPkt, fs, -fineCFO);


fprintf('Coarse CFO=%.1f Hz, Fine CFO=%.1f Hz\n', coarseCFO, fineCFO);

fprintf('pktStart=%d, symOffset=%d\n', pktStart, symOffset);

% Demod L-LTF and estimate channel
ofdmInfo = wlanNonHTOFDMInfo('NonHT-Data', cfgNonHT);
%disp(ofdmInfo)
% OFDM info for NonHT-Data

Nfft = ofdmInfo.FFTLength;
Ncp  = ofdmInfo.CPLength;          % your struct uses CPLength field name
activeFFT = ofdmInfo.ActiveFFTIndices;   % 52 indices into FFT bins (MATLAB 1..Nfft)
pilotIdx52 = ofdmInfo.PilotIndices;      % 4 indices into the 52-tone active vector
dataIdx52  = ofdmInfo.DataIndices;       % 48 indices into the 52-tone active vector

symLen = Nfft + Ncp;

rxData = yPkt(ind.NonHTData(1):ind.NonHTData(2), :);

fprintf('rxData length: %d samples\n', size(rxData,1));

numSym = floor(size(rxData,1) / symLen);

% Extract L-LTF from received packet
rxLLTF = yPkt(ind.LLTF(1):ind.LLTF(2), :);

% Demod + channel estimate + noise estimate
lltfDemod = wlanLLTFDemodulate(rxLLTF, cfgNonHT);
chanEst = wlanLLTFChannelEstimate(lltfDemod, cfgNonHT);
noiseEst = wlanLLTFNoiseEstimate(lltfDemod);

% Pull channel estimate on active tones (assume 1 RX antenna)
H = squeeze(chanEst(:,1,1));   % should be 52x1
assert(numel(H)==numel(activeFFT), 'chanEst length does not match NumTones(52).');

rxPilotsEq = zeros(numel(pilotIdx52), numSym);
rxDataEq   = zeros(numel(dataIdx52),  numSym);

for s = 1:numSym
    seg = rxData((s-1)*symLen + (1:symLen), 1);   % one OFDM symbol time-domain (1 Rx)
    seg = seg(Ncp+1:end);                         % remove CP
    Y = fft(seg, Nfft);                           % FFT bins 1..Nfft

    % Select the 52 active tones in WLAN ordering
    Yact = Y(activeFFT);

    % Equalize
    Zeq = Yact ./ (H + 1e-12);

    % Extract pilots/data in the 52-tone domain
    rxPilotsEq(:,s) = Zeq(pilotIdx52);
    rxDataEq(:,s)   = Zeq(dataIdx52);
end

% Plot equalized pilots (all symbols)
figure;
plot(real(rxPilotsEq(:)), imag(rxPilotsEq(:)), '.'); grid on;
title('Equalized pilots (NonHT-Data) across symbols'); xlabel('I'); ylabel('Q');

% Field indices for Non-HT packet
ind = wlanFieldIndices(cfgNonHT);

% Extract Non-HT DATA field from received packet
if size(yPkt,1) < ind.NonHTData(2)
    error('yPkt too short after timing correction; increase SamplesPerFrame or margin.');
end

% OFDM demodulate/equalize DATA symbols
%nonHTDataDemod = wlanNonHTDataDemodulate(rxData, cfgNonHT, chanEst);
% Recover PSDU and get equalized data symbols
[rxPSDU, eqSym, rxState] = wlanNonHTDataRecover(rxData, chanEst, noiseEst, cfgNonHT);

fprintf('eqSym size = %dx%d, rxDataEq size = %dx%d\n', size(eqSym,1), size(eqSym,2), size(rxDataEq,1), size(rxDataEq,2));
Ncmp = min(size(eqSym,2), size(rxDataEq,2));
err = rxDataEq(:,1:Ncmp) - eqSym(:,1:Ncmp);
fprintf('Mean abs diff between rxDataEq and eqSym: %.3e\n', mean(abs(err(:))));

z = eqSym(:);
% z = z(1:min(5000,numel(z)));
% figure; plot(real(z), imag(z), '.'); grid on;
% title('Equalized constellation (subset)');
% xlabel('I'); ylabel('Q');
% Random subsample for plotting clarity
Nplot = min(30000, numel(z));
idx = randperm(numel(z), Nplot);
zplot = z(idx);

figure; plot(real(zplot), imag(zplot), '.'); grid on;
title('Equalized DATA constellation (random subset)');
xlabel('I'); ylabel('Q');

%figure; histogram(real(zplot), 100); grid on; title('Histogram of I');
%figure; histogram(imag(zplot), 100); grid on; title('Histogram of Q');

figure; histogram(real(eqSym(:)), 200); grid on;
title('Histogram of real(eqSym) (MCS0 BPSK)');

z = eqSym(:);

% BPSK slicer (MCS 0)
% shat = sign(real(z));
% shat(shat==0) = 1;      % avoid zeros
% shat = complex(shat, 0);
% 
% % Best-fit complex scalar gain
% g = (shat' * z) / (shat' * shat);
% 
% % RMS EVM
% evm_rms = sqrt(mean(abs(z - g*shat).^2) / mean(abs(g*shat).^2));
% fprintf('EVM_rms = %.3f (%.2f dB)\n', evm_rms, 20*log10(evm_rms));
z = eqSym(:);

% QPSK slicer
shat = sign(real(z)) + 1j*sign(imag(z));
shat = shat / sqrt(2);

% Fit complex gain
g = (shat' * z) / (shat' * shat);

evm = sqrt(mean(abs(z - g*shat).^2) / mean(abs(g*shat).^2));
fprintf('QPSK EVM_rms = %.3f (%.2f dB)\n', evm, 20*log10(evm));

% Trimmed EVM to reduce influence of outliers (optional)
e = abs(z - g*shat).^2;
e_sorted = sort(e);
keep = round(0.95*numel(e));
evm_trim = sqrt(mean(e_sorted(1:keep)) / mean(abs(g*shat).^2));
fprintf('EVM_trim(95%%) = %.3f (%.2f dB)\n', evm_trim, 20*log10(evm_trim));


% Recover non-HT Data field

fprintf('Noise estimate: %.3e\n', noiseEst);

pktLen = ind.NonHTData(2);
if size(yPkt,1) < pktLen
    error('Captured packet too short: need %d samples after pktStart, got %d. Increase SamplesPerFrame.', ...
        pktLen, size(yPkt,1));
end

% Constellation sanity (just show a subset)
% figure;
% plot(real(nonHTDataDemod(:)), imag(nonHTDataDemod(:)), '.');
% grid on; title('Equalized non-HT DATA constellation (WLAN toolbox)');
% xlabel('I'); ylabel('Q');

fprintf('NonHTDataRecover: num bits recovered = %d\n', numel(rxPSDU));

%% ------------ Your authentication metrics placeholder ------------
% At this stage, you can compute perm/IM metrics on pilots or selected tones.
% For A2-1, the "pilot sign coding" is deterministic; receiver can predict it given seed.
% Next step: extract pilots per OFDM symbol and compute your correlator/IM stats.
disp('Receiver chain complete. Next: extract pilots per symbol and compute metrics.');
% psdu and rxPSDU should be column vectors of bits (0/1)
fprintf('TX psdu bits: %d | RX psdu bits: %d\n', numel(psdu), numel(rxPSDU));

numErr = sum(rxPSDU ~= psdu);
ber = numErr / numel(psdu);
L = min(numel(psdu), numel(rxPSDU));
numErr = sum(rxPSDU(1:L) ~= psdu(1:L));
ber = numErr / L;
fprintf('Bit errors: %d/%d (BER=%.3e)\n', numErr, L, ber);
%fprintf('Bit errors: %d / %d (BER=%.3e)\n', numErr, numel(psdu), ber);


%% ==========================================================
function y = applyPilotSignCodingNonHT(waveform, cfgNonHT, seed)
% Placeholder "pilot coding" hook.
% The easiest safe A2-1: multiply the *entire time-domain waveform* by a slow phase
% is NOT correct. We need to modify pilots in frequency domain per OFDM symbol,
% but WLAN Toolbox abstracts that.
%
% So for now, this function does NOTHING.
%
% Tomorrow: we will implement pilot modification by regenerating the waveform
% with custom OFDM symbol mapping (or by patching the DATA OFDM symbols).
y = waveform;

% Keep deterministic seed visible so we don't forget
rng(double(seed));
end