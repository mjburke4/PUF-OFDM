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
% waveform = wlanWaveformGenerator(psdu, cfgNonHT, ...
%     'IdleTime', idleTime);
waveform_uncoded = wlanWaveformGenerator(psdu, cfgNonHT, ...
    'IdleTime', idleTime);

% Extract the nominal TX payload pilots BEFORE coding
txPilotsBase = extractPayloadPilotsFromWaveform(waveform_uncoded, cfgNonHT);
fprintf('txPilotsBase size: %dx%d\n', size(txPilotsBase,1), size(txPilotsBase,2));
waveform = waveform_uncoded;

% Optional: apply pilot sign coding to DATA symbols (A2-1)
% if enablePilotCoding
%     waveform = applyPilotSignCodingNonHT(waveform, cfgNonHT, seed);
% end
% ---- Base challenge + frame nonce for this packet ----
B = 64;
%% ------------ PUF node setup ------------
U = 1;              % one enrolled transmitter for now
u0 = 1;             % legitimate node index
m = 64;             % challenge length
B = 64;             % number of PUF bits used for seed derivation

nodes = init_node_population(U, m, 'puf_bits', 128, 'seed', 7);

node = nodes(u0);          % or however your enrolled node is defined
node_claimed = node;
%C_seed = randi([0 1], 1, nodes(u0).m);
m = 64;
%C_seed = randi([0 1],1,m);
C_seed = randi([0 1], 1, node.m);
nonce_frame = randi([0 1], 1, 32);

if enablePilotCoding
    %waveform = applyPilotSignCodingNonHT(waveform, cfgNonHT, nodes(u0), C_seed, nonce_frame, B);
   %waveform = applyPilotSignCodingNonHT(waveform, cfgNonHT, C_seed, nonce_frame, B);
    waveform = applyPilotSignCodingNonHT(waveform, cfgNonHT, node, C_seed, nonce_frame, B);
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

% ---- Extract equalized pilots and data per OFDM symbol from rxData ----
ofdmInfo = wlanNonHTOFDMInfo('NonHT-Data', cfgNonHT);

Nfft = ofdmInfo.FFTLength;
Ncp  = ofdmInfo.CPLength;
activeFFT = ofdmInfo.ActiveFFTIndices;   % 52 active FFT bins
pilotIdx52 = ofdmInfo.PilotIndices;      % 4 pilot indices within 52 active tones
dataIdx52  = ofdmInfo.DataIndices;       % 48 data indices within 52 active tones

symLen = Nfft + Ncp;
numSym = floor(size(rxData,1)/symLen);

H = squeeze(chanEst(:,1,1));  % 52x1 for single RX antenna

rxPilotsEq = zeros(numel(pilotIdx52), numSym);
rxDataEq   = zeros(numel(dataIdx52), numSym);

for ell = 1:numSym
    seg = rxData((ell-1)*symLen + (1:symLen), 1);
    seg = seg(Ncp+1:end);         % remove CP
    Y = fft(seg, Nfft);           % 64 FFT bins

    % Select 52 active tones
    Yact = Y(activeFFT);

    % Equalize
    Zeq = Yact ./ (H + 1e-12);

    % Split pilots and data
    rxPilotsEq(:,ell) = Zeq(pilotIdx52);
    rxDataEq(:,ell)   = Zeq(dataIdx52);
end

fprintf('Extracted pilots: %dx%d\n', size(rxPilotsEq,1), size(rxPilotsEq,2));
fprintf('Extracted data  : %dx%d\n', size(rxDataEq,1), size(rxDataEq,2));

figure;
plot(real(rxPilotsEq(:)), imag(rxPilotsEq(:)), '.');
grid on;
title('Equalized received pilots across OFDM symbols');
xlabel('I'); ylabel('Q');

figure;
for k = 1:4
    subplot(4,1,k);
    stem(real(rxPilotsEq(k,:)), 'filled');
    grid on;
    ylabel(sprintf('P%d',k));
end
sgtitle('Real part of equalized payload pilots vs OFDM symbol index');
xlabel('OFDM symbol index');

% % ---- Phase-coded pilot correlator metric (matches TX model) ----
% 
% Tp_H1 = zeros(1,numSym);
% Tp_H0 = zeros(1,numSym);
% phiH1_deg = zeros(1,numSym);
% phiH0_deg = zeros(1,numSym);
% 
% % Baseline pilot reference.
% % Start with all-ones since we are measuring relative phase-coded coherence.
% p_base = ones(numel(pilotIdx52),1);
% 
% for ell = 1:numSym
%     z = rxPilotsEq(:,ell);
% 
%     % ----- H1: correct seed rule -----
%     nonce_bits_H1 = nonce_with_symbol(nonce_frame, ell, 16);
%     C_eff_H1      = mix_challenge_nonce(C_seed, nonce_bits_H1, length(C_seed));
% 
%     % Must match TX-side seed rule exactly
%     seed_H1 = uint32(mod(sum(double(C_eff_H1)) + 7919*ell + 104729, 2^32-1));
% 
%     phi1_deg = phase_from_seed(seed_H1);
%     phi1 = deg2rad(phi1_deg);
% 
%     p1 = p_base * exp(1j*phi1);
% 
%     % ----- H0: wrong/stale seed rule -----
%     % Simple wrong hypothesis: perturb the frame nonce before reconstruction
%     nonce_frame_H0 = [1 - nonce_frame(1), nonce_frame(2:end)];
%     nonce_bits_H0  = nonce_with_symbol(nonce_frame_H0, ell, 16);
%     C_eff_H0       = mix_challenge_nonce(C_seed, nonce_bits_H0, length(C_seed));
% 
%     seed_H0 = uint32(mod(sum(double(C_eff_H0)) + 7919*ell + 104729, 2^32-1));
% 
%     phi0_deg = phase_from_seed(seed_H0);
%     phi0 = deg2rad(phi0_deg);
% 
%     p0 = p_base * exp(1j*phi0);
% 
%     % Coherence-sensitive correlator
%     Tp_H1(ell) = real((p1' * z) / (norm(p1)^2 + 1e-12));
%     Tp_H0(ell) = real((p0' * z) / (norm(p0)^2 + 1e-12));
% 
%     phiH1_deg(ell) = phi1_deg;
%     phiH0_deg(ell) = phi0_deg;
% end
% 
% figure;
% histogram(Tp_H1, 40, 'Normalization','pdf');
% hold on;
% histogram(Tp_H0, 40, 'Normalization','pdf');
% grid on;
% legend('H1 (correct phase hypothesis)','H0 (wrong phase hypothesis)');
% title('Phase-coded pilot correlator statistic: H1 vs H0');
% xlabel('T_p'); ylabel('PDF');
% 
% fprintf('Mean Tp_H1 = %.3f\n', mean(Tp_H1));
% fprintf('Mean Tp_H0 = %.3f\n', mean(Tp_H0));
% fprintf('Std  Tp_H1 = %.3f\n', std(Tp_H1));
% fprintf('Std  Tp_H0 = %.3f\n', std(Tp_H0));
% 
% fprintf('First 10 H1 phases (deg): ');
% fprintf('%+g ', phiH1_deg(1:min(10,numSym)));
% fprintf('\n');
% 
% fprintf('First 10 H0 phases (deg): ');
% fprintf('%+g ', phiH0_deg(1:min(10,numSym)));
% fprintf('\n');

% ---- Phase-hypothesis detector for common pilot rotation ----

% Tp_H1 = zeros(1,numSym);
% Tp_H0 = zeros(1,numSym);
% angErr_H1_deg = zeros(1,numSym);
% angErr_H0_deg = zeros(1,numSym);
% 
% phiH1_deg = zeros(1,numSym);
% phiH0_deg = zeros(1,numSym);
% 
% for ell = 1:numSym
%     z = rxPilotsEq(:,ell);
%     p_nom = txPilotsBase(:,ell);
% 
%     % Collapse received pilots onto nominal TX pilot direction
%     u = (p_nom' * z) / (norm(p_nom)^2 + 1e-12);
% 
%     % ----- H1 -----
%     nonce_bits_H1 = nonce_with_symbol(nonce_frame, ell, 16);
%     C_eff_H1      = mix_challenge_nonce(C_seed, nonce_bits_H1, length(C_seed));
%     seed_H1       = uint32(mod(sum(double(C_eff_H1)) + 7919*ell + 104729, 2^32-1));
%     phi1_deg      = phase_from_seed(seed_H1);
%     phi1          = deg2rad(phi1_deg);
% 
%     % ----- H0 -----
%     nonce_frame_H0 = [1 - nonce_frame(1), nonce_frame(2:end)];
%     nonce_bits_H0  = nonce_with_symbol(nonce_frame_H0, ell, 16);
%     C_eff_H0       = mix_challenge_nonce(C_seed, nonce_bits_H0, length(C_seed));
%     seed_H0        = uint32(mod(sum(double(C_eff_H0)) + 7919*ell + 104729, 2^32-1));
%     phi0_deg       = phase_from_seed(seed_H0);
%     phi0           = deg2rad(phi0_deg);
% 
%     % Residual angle between measured pilot phasor and each hypothesis
%     e1 = angle(u * exp(-1j*phi1));
%     e0 = angle(u * exp(-1j*phi0));
% 
%     % Scores: larger is better, max = 1
%     Tp_H1(ell) = cos(e1);
%     Tp_H0(ell) = cos(e0);
% 
%     angErr_H1_deg(ell) = rad2deg(e1);
%     angErr_H0_deg(ell) = rad2deg(e0);
% 
%     phiH1_deg(ell) = phi1_deg;
%     phiH0_deg(ell) = phi0_deg;
% end
% 
% figure;
% histogram(Tp_H1, 40, 'Normalization','pdf');
% hold on;
% histogram(Tp_H0, 40, 'Normalization','pdf');
% grid on;
% legend('H1','H0');
% title('Phase-consistency score: H1 vs H0');
% xlabel('cos(angle error)'); ylabel('PDF');
% 
% fprintf('Mean Tp_H1 = %.3f\n', mean(Tp_H1));
% fprintf('Mean Tp_H0 = %.3f\n', mean(Tp_H0));
% fprintf('Std  Tp_H1 = %.3f\n', std(Tp_H1));
% fprintf('Std  Tp_H0 = %.3f\n', std(Tp_H0));
% 
% fprintf('Mean |angle err| H1 = %.2f deg\n', mean(abs(angErr_H1_deg)));
% fprintf('Mean |angle err| H0 = %.2f deg\n', mean(abs(angErr_H0_deg)));
% 
% fprintf('First 10 H1 phases (deg): ');
% fprintf('%+g ', phiH1_deg(1:min(10,numSym)));
% fprintf('\n');
% 
% fprintf('First 10 H0 phases (deg): ');
% fprintf('%+g ', phiH0_deg(1:min(10,numSym)));
% fprintf('\n');
% 
% figure;
% plot(rad2deg(angle((sum(conj(txPilotsBase).*rxPilotsEq,1)) ./ ...
%     (sum(abs(txPilotsBase).^2,1) + 1e-12))), 'o-');
% grid on;
% title('Estimated common pilot rotation per OFDM symbol');
% xlabel('OFDM symbol index');
% ylabel('Estimated phase (deg)');

% ---- Baseline / coded pilot metric ----

% Standard reference pilot vector (placeholder baseline)
% Start simple: assume all-ones pilot reference.

% p_base = ones(4,1);
% 
% % Deterministic sign sequence for H1 and H0
% seed_H1 = 12345;
% seed_H0 = 54321;
% 
% S_H1 = pilotSignSequence(numSym, seed_H1);   % correct seed
% S_H0 = pilotSignSequence(numSym, seed_H0);   % wrong seed
% 
% Tp_base = zeros(1,numSym);
% Tp_H1   = zeros(1,numSym);
% Tp_H0   = zeros(1,numSym);
% 
% for ell = 1:numSym
%     z = rxPilotsEq(:,ell);
% 
%     p1 = p_base .* S_H1(:,ell);   % expected pilot signs under H1
%     p0 = p_base .* S_H0(:,ell);   % expected pilot signs under H0
% 
%     Tp_base(ell) = real((p_base' * z) / (p_base' * p_base));
%     Tp_H1(ell)   = real((p1'    * z) / (p1'    * p1));
%     Tp_H0(ell)   = real((p0'    * z) / (p0'    * p0));
% end

% figure;
% histogram(Tp_H1, 40, 'Normalization','pdf');
% hold on;
% histogram(Tp_H0, 40, 'Normalization','pdf');
% grid on;
% legend('H1 (correct seed)','H0 (wrong seed)');
% title('Pilot correlator statistic: H1 vs H0');
% xlabel('T_p'); ylabel('PDF');
% 
% fprintf('Mean Tp_H1 = %.3f\n', mean(Tp_H1));
% fprintf('Mean Tp_H0 = %.3f\n', mean(Tp_H0));
% fprintf('Std  Tp_H1 = %.3f\n', std(Tp_H1));
% fprintf('Std  Tp_H0 = %.3f\n', std(Tp_H0));
% ---- Differential phase detector across OFDM symbols ----

u_meas = zeros(1,numSym);
phiH1_deg = zeros(1,numSym);
phiH0_deg = zeros(1,numSym);

for ell = 1:numSym
    z = rxPilotsEq(:,ell);
    p_nom = txPilotsBase(:,ell);

    % Measured common pilot phasor relative to nominal TX pilot template
    u_meas(ell) = (p_nom' * z) / (norm(p_nom)^2 + 1e-12);

    % ----- H1 phase sequence -----
    % nonce_bits_H1 = nonce_with_symbol(nonce_frame, ell, 16);
    % C_eff_H1      = mix_challenge_nonce(C_seed, nonce_bits_H1, length(C_seed));
    % seed_H1       = uint32(mod(sum(double(C_eff_H1)) + 7919*ell + 104729, 2^32-1));
    % phiH1_deg(ell)= phase_from_seed(seed_H1);
    nonce_bits_H1 = nonce_with_symbol(nonce_frame, ell, 16);
    C_eff_H1      = mix_challenge_nonce(C_seed, nonce_bits_H1, length(C_seed));
    
    R_H1          = get_puf_bits(node_claimed, C_eff_H1, B);
    seed_H1       = hash_bits_u32(R_H1);
    
    phiH1_deg(ell)= phase_from_seed(seed_H1);

    % ----- H0 phase sequence -----
    % nonce_frame_H0 = [1 - nonce_frame(1), nonce_frame(2:end)];
    % nonce_bits_H0  = nonce_with_symbol(nonce_frame_H0, ell, 16);
    % C_eff_H0       = mix_challenge_nonce(C_seed, nonce_bits_H0, length(C_seed));
    % seed_H0        = uint32(mod(sum(double(C_eff_H0)) + 7919*ell + 104729, 2^32-1));
    % phiH0_deg(ell) = phase_from_seed(seed_H0);
    nonce_frame_H0 = [1 - nonce_frame(1), nonce_frame(2:end)];
    nonce_bits_H0  = nonce_with_symbol(nonce_frame_H0, ell, 16);
    C_eff_H0       = mix_challenge_nonce(C_seed, nonce_bits_H0, length(C_seed));
    
    R_H0           = get_puf_bits(node_claimed, C_eff_H0, B);
    seed_H0        = hash_bits_u32(R_H0);
    
    phiH0_deg(ell) = phase_from_seed(seed_H0);
end

% Differential scores begin at ell = 2
Td_H1 = zeros(1,numSym-1);
Td_H0 = zeros(1,numSym-1);
dErr_H1_deg = zeros(1,numSym-1);
dErr_H0_deg = zeros(1,numSym-1);

for ell = 2:numSym
    % Measured differential pilot phase
    d_meas = u_meas(ell) * conj(u_meas(ell-1));
    dphi_meas = angle(d_meas);

    % Hypothesized differential phases
    dphi_H1 = deg2rad(phiH1_deg(ell) - phiH1_deg(ell-1));
    dphi_H0 = deg2rad(phiH0_deg(ell) - phiH0_deg(ell-1));

    % Wrap residual errors to [-pi, pi]
    e1 = angle(exp(1j*(dphi_meas - dphi_H1)));
    e0 = angle(exp(1j*(dphi_meas - dphi_H0)));

    % Cosine score: +1 is best
    Td_H1(ell-1) = cos(e1);
    Td_H0(ell-1) = cos(e0);

    dErr_H1_deg(ell-1) = rad2deg(e1);
    dErr_H0_deg(ell-1) = rad2deg(e0);
end

figure;
histogram(Td_H1, 40, 'Normalization','pdf');
hold on;
histogram(Td_H0, 40, 'Normalization','pdf');
grid on;
legend('H1','H0');
title('Differential phase score: H1 vs H0');
xlabel('cos(differential phase error)'); ylabel('PDF');

fprintf('Mean Td_H1 = %.3f\n', mean(Td_H1));
fprintf('Mean Td_H0 = %.3f\n', mean(Td_H0));
fprintf('Std  Td_H1 = %.3f\n', std(Td_H1));
fprintf('Std  Td_H0 = %.3f\n', std(Td_H0));

fprintf('Mean |diff angle err| H1 = %.2f deg\n', mean(abs(dErr_H1_deg)));
fprintf('Mean |diff angle err| H0 = %.2f deg\n', mean(abs(dErr_H0_deg)));


fprintf('First 10 H1 phases (deg): ');
fprintf('%+g ', phiH1_deg(1:min(10,numSym)));
fprintf('\n');

fprintf('First 10 H0 phases (deg): ');
fprintf('%+g ', phiH0_deg(1:min(10,numSym)));
fprintf('\n');

%% ---- DPPA authentication decision layer ----
% Differential detector gives numSym-1 scores because it uses adjacent symbols

Ldiff = numel(Td_H1);

% Start with a hand-tuned threshold and vote ratio
tau_d = 0.95;      % per-differential-symbol threshold
alpha = 0.80;      % frame acceptance ratio

% Per-symbol decisions
b_H1 = (Td_H1 >= tau_d);
b_H0 = (Td_H0 >= tau_d);

% Vote counts
V_H1 = sum(b_H1);
V_H0 = sum(b_H0);

% Frame decisions
accept_H1 = (V_H1 >= ceil(alpha * Ldiff));
accept_H0 = (V_H0 >= ceil(alpha * Ldiff));

fprintf('\n=== DPPA Frame Authentication Decision ===\n');
fprintf('Ldiff = %d, tau_d = %.3f, alpha = %.2f\n', Ldiff, tau_d, alpha);
fprintf('H1 votes = %d / %d  --> accept_H1 = %d\n', V_H1, Ldiff, accept_H1);
fprintf('H0 votes = %d / %d  --> accept_H0 = %d\n', V_H0, Ldiff, accept_H0);

figure;
stem(1:Ldiff, Td_H1, 'filled'); hold on;
yline(tau_d, 'r--', 'LineWidth', 1.5);
grid on;
title('DPPA differential scores under H1 with threshold');
xlabel('Differential symbol index');
ylabel('T_d');

%% ---- Quick threshold sweep on current packet ----
tau_vec = 0.70:0.01:0.999;

passFrac_H1 = zeros(size(tau_vec));
passFrac_H0 = zeros(size(tau_vec));

for k = 1:numel(tau_vec)
    tau_k = tau_vec(k);
    passFrac_H1(k) = mean(Td_H1 >= tau_k);
    passFrac_H0(k) = mean(Td_H0 >= tau_k);
end

figure;
plot(tau_vec, passFrac_H1, 'LineWidth', 1.5); hold on;
plot(tau_vec, passFrac_H0, 'LineWidth', 1.5);
grid on;
legend('H1 pass fraction','H0 pass fraction','Location','best');
title('DPPA symbol-pass fraction vs threshold');
xlabel('\tau_d');
ylabel('Fraction of differential symbols passing');

figure;
stem(1:Ldiff, Td_H0, 'filled'); hold on;
yline(tau_d, 'r--', 'LineWidth', 1.5);
grid on;
title('DPPA differential scores under H0 with threshold');
xlabel('Differential symbol index');
ylabel('T_d');

figure;
plot(rad2deg(angle(u_meas)), 'o-');
grid on;
title('Measured pilot phase per OFDM symbol');
xlabel('OFDM symbol index');
ylabel('Phase (deg)');

figure;
plot(rad2deg(angle(u_meas(2:end).*conj(u_meas(1:end-1)))), 'o-');
grid on;
title('Measured differential pilot phase');
xlabel('Differential symbol index');
ylabel('Phase difference (deg)');

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
% function y = applyPilotSignCodingNonHT(waveform, cfgNonHT, seed)
% % Placeholder "pilot coding" hook.
% % The easiest safe A2-1: multiply the *entire time-domain waveform* by a slow phase
% % is NOT correct. We need to modify pilots in frequency domain per OFDM symbol,
% % but WLAN Toolbox abstracts that.
% %
% % So for now, this function does NOTHING.
% %
% % Tomorrow: we will implement pilot modification by regenerating the waveform
% % with custom OFDM symbol mapping (or by patching the DATA OFDM symbols).
% y = waveform;
% 
% % Keep deterministic seed visible so we don't forget
% rng(double(seed));
% end
%function y = applyPilotSignCodingNonHT(waveform, cfgNonHT, C_seed, nonce_frame, B)% applyPilotSignCodingNonHT
function y = applyPilotSignCodingNonHT(waveform, cfgNonHT, node, C_seed, nonce_frame, B)
% Apply deterministic per-symbol pilot sign coding to the NonHT-Data field
% only, leaving preamble/training/SIG untouched.
%
% Inputs:
%   waveform     : complex baseband WLAN packet from wlanWaveformGenerator
%   cfgNonHT     : wlanNonHTConfig
%   node         : enrolled node struct (must support node.m and get_puf_bits)
%   C_seed       : base challenge bits, row/col OK
%   nonce_frame  : frame nonce bits, row/col OK
%   B            : number of PUF bits used to derive seed
%
% Output:
%   y            : patched waveform, same size as waveform

    y = waveform;

    % Field/sample indices within ONE WLAN packet
    ind = wlanFieldIndices(cfgNonHT);

    % OFDM info for NonHT-Data
    ofdmInfo   = wlanNonHTOFDMInfo('NonHT-Data', cfgNonHT);
    Nfft       = ofdmInfo.FFTLength;
    Ncp        = ofdmInfo.CPLength;
    symLen     = Nfft + Ncp;
    activeFFT  = ofdmInfo.ActiveFFTIndices;   % 52 active bins in FFT domain
    pilotIdx52 = ofdmInfo.PilotIndices;       % 4 pilot positions within 52-tone active set

    dataStart = ind.NonHTData(1);
    dataEnd   = ind.NonHTData(2);

    rxData = y(dataStart:dataEnd, 1);

    if mod(numel(rxData), symLen) ~= 0
        error('NonHT-Data length (%d) is not an integer multiple of OFDM symbol length (%d).', ...
            numel(rxData), symLen);
    end

    numSym = numel(rxData) / symLen;

    for ell = 1:numSym
        % ----- Extract one payload OFDM symbol -----
        idx = (ell-1)*symLen + (1:symLen);
        xcp = rxData(idx);

        cp  = xcp(1:Ncp);
        x   = xcp(Ncp+1:end);

        % FFT to 64-bin domain
        X64 = fft(x, Nfft);

        % Pull 52 active tones
        X52 = X64(activeFFT);

        % ----- Derive per-symbol sign pattern from your PUF pipeline -----
        % nonce_bits = nonce_with_symbol(nonce_frame, ell, 16);
        % C_eff = mix_challenge_nonce(C_seed, nonce_bits, length(C_seed));
        % 
        % R_bits     = get_puf_bits(node, C_eff, B);
        % seed_u32   = hash_bits_u32(R_bits);
        % 
        % s = pilot_signs_from_seed(seed_u32, numel(pilotIdx52));   % 4x1, entries ±1
        % nonce_bits = nonce_with_symbol(nonce_frame, ell, 16);
        % C_eff      = mix_challenge_nonce(C_seed, nonce_bits, length(C_seed));
        % 
        % % Temporary deterministic non-PUF seed for OTA bring-up
        % seed_u32 = uint32(mod(sum(double(C_eff)) + 7919*ell + 104729, 2^32-1));
        % 
        % phi_deg = phase_from_seed(seed_u32);
        nonce_bits = nonce_with_symbol(nonce_frame, ell, 16);
        C_eff      = mix_challenge_nonce(C_seed, nonce_bits, length(C_seed));
        
        % Real PUF-derived seed
        R_bits   = get_puf_bits(node, C_eff, B);
        %seed_u32 = hash_bits_u32(R_bits);
        R_bits   = get_puf_bits(node, C_eff, B);
        seed_u32 = hash_bits_u32(R_bits);
        
        phi_deg = phase_from_seed(seed_u32);
        phi = deg2rad(phi_deg);

        rot = exp(1j*phi);

        if ell <= 3
            fprintf('TX symbol %d seed=%u phi=%+.1f deg\n', ...
                ell, seed_u32, phi_deg);
        end

        X52(pilotIdx52) = X52(pilotIdx52) .* rot;

        % ----- Write back into full FFT vector and reconstruct symbol -----
        X64(activeFFT) = X52;
        x_mod = ifft(X64, Nfft);

        % Restore CP
        xcp_mod = [x_mod(end-Ncp+1:end); x_mod];

        % Overwrite this symbol in payload region
        rxData(idx) = xcp_mod;
        % if ell <= 3
        %     fprintf('TX symbol %d seed=%u signs=[%+d %+d %+d %+d]\n', ...
        %         ell, seed_u32, s(1), s(2), s(3), s(4));
        % end
    end

    % Write patched payload back into output waveform
    y(dataStart:dataEnd, 1) = rxData;
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

function phi_deg = phase_from_seed(seed_u32)
% Must match TX-side phase mapping exactly
    %phi_choices_deg = [-24, -12, +12, +24];
    phi_choices_deg = [-8, +8];
    idx_phi = 1 + mod(double(seed_u32), numel(phi_choices_deg));
    phi_deg = phi_choices_deg(idx_phi);
end

function txPilots = extractPayloadPilotsFromWaveform(waveform, cfgNonHT)
% Extract payload pilot vectors from an uncoded WLAN waveform

    ind = wlanFieldIndices(cfgNonHT);

    ofdmInfo   = wlanNonHTOFDMInfo('NonHT-Data', cfgNonHT);
    Nfft       = ofdmInfo.FFTLength;
    Ncp        = ofdmInfo.CPLength;
    symLen     = Nfft + Ncp;
    activeFFT  = ofdmInfo.ActiveFFTIndices;
    pilotIdx52 = ofdmInfo.PilotIndices;

    xData = waveform(ind.NonHTData(1):ind.NonHTData(2), 1);

    if mod(numel(xData), symLen) ~= 0
        error('Payload length is not an integer multiple of OFDM symbol length.');
    end

    numSym = numel(xData) / symLen;
    txPilots = zeros(numel(pilotIdx52), numSym);

    for ell = 1:numSym
        idx = (ell-1)*symLen + (1:symLen);
        xcp = xData(idx);
        x   = xcp(Ncp+1:end);

        X64 = fft(x, Nfft);
        X52 = X64(activeFFT);

        txPilots(:,ell) = X52(pilotIdx52);
    end
end
