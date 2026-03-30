%% Pluto OFDM synthetic preamble bring-up (2 Plutos)
clear; clc;

%% -------- User settings --------
fc = 2.462e9;          % center freq (try ch 11)
fs = 1e6;              % sample rate (Hz)
txGain = -10;          % dB (start low)
rxGain = 30;           % dB

Nfft = 64;
Ncp  = 16;
Nsym = 10;             % OFDM payload symbols per frame (keep small for now)

Np = 8;                % number of pilots (your model)
pilotIdx = round(linspace(2, Nfft-1, Np));  % simple spread
pilotIdx = unique(pilotIdx); pilotIdx = pilotIdx(1:Np);

% Fixed, distinct pilot values (unit magnitude)
pilots = exp(1j*2*pi*(0:Np-1)/Np).';

% Preamble: repeat a short block (timing-friendly)
Lshort = 16;
shortBlock = exp(1j*2*pi*(0:Lshort-1)/Lshort).'; % constant-modulus
Nrep = 10;                                       % repetition count
preamble = repmat(shortBlock, Nrep, 1);

% Frame spacing
guard = zeros(200,1);     % quiet gap for easier sync
%% --------------------------------

%% Build OFDM payload (time-domain)
payload_td = zeros((Nfft+Ncp)*Nsym, 1);

for s = 1:Nsym
    X = zeros(Nfft,1);

    % Put pilots
    X(pilotIdx) = pilots;

    % Fill other bins with QPSK
    dataIdx = setdiff((1:Nfft).', pilotIdx(:));
    bits = randi([0 1], 2*numel(dataIdx), 1);
    qpsk = (2*bits(1:2:end)-1) + 1j*(2*bits(2:2:end)-1);
    qpsk = qpsk / sqrt(2);

    X(dataIdx) = qpsk;

    x = ifft(X, Nfft);
    x_cp = [x(end-Ncp+1:end); x];

    payload_td((s-1)*(Nfft+Ncp)+1 : s*(Nfft+Ncp)) = x_cp;
end

% Full frame: guard + preamble + guard + payload + guard
frame = [guard; preamble; guard; payload_td; guard];

% Normalize to avoid clipping
frame = 0.6 * frame / max(abs(frame)+1e-12);

%% Pluto TX (usb:0)
tx = sdrtx('Pluto', ...
    'RadioID','usb:0', ...
    'CenterFrequency',fc, ...
    'BasebandSampleRate',fs, ...
    'Gain',txGain);

%% Pluto RX (usb:1)
rx = sdrrx('Pluto', ...
    'RadioID','usb:1', ...
    'CenterFrequency',fc, ...
    'BasebandSampleRate',fs, ...
    'SamplesPerFrame', 200000, ...
    'GainSource','Manual', ...
    'Gain',rxGain, ...
    'OutputDataType','double');

%% Transmit repeat + receive one capture
disp('Starting transmitRepeat...');
transmitRepeat(tx, frame);

pause(0.2); % let it settle

disp('Capturing...');
y = rx();

% Stop
release(tx);
release(rx);

%% Sync: correlate with preamble
r = abs(conv(y, flipud(conj(preamble)), 'valid'));
[~, idx] = max(r);
start = idx;  % estimated start of preamble in y

fprintf('Estimated preamble start index: %d\n', start);

% Extract one OFDM symbol after preamble+guard
preambleLen = numel(preamble);
guardLen = numel(guard);

payloadStart = start + preambleLen + guardLen; % where payload should begin
sym1 = y(payloadStart : payloadStart + (Nfft+Ncp) - 1);

% Remove CP and FFT
sym1_nocp = sym1(Ncp+1:end);
Y1 = fft(sym1_nocp, Nfft);

% Plot constellation of data bins
dataIdx = setdiff((1:Nfft).', pilotIdx(:));
figure; plot(real(Y1(dataIdx)), imag(Y1(dataIdx)), '.'); grid on;
title('Raw constellation (no EQ)'); xlabel('I'); ylabel('Q');

% Plot sync correlation peak
figure; plot(r); grid on; title('Preamble correlation magnitude'); xlabel('Sample'); ylabel('|corr|');

disp('Done.');