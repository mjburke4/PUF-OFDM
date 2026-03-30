% ===== USER SETTINGS =====
fc = 2.45e9;          % center frequency
fs = 1e6;             % sample rate
txGain = -10;         % start low for OTA
rxGain = 30;          % manual gain to start
toneHz = 50e3;        % tone frequency
N = 200000;           % samples for transmit/receive
% =========================

% --- TX Pluto (usb:0) ---
tx = sdrtx('Pluto', ...
    'RadioID','usb:0', ...
    'CenterFrequency',fc, ...
    'BasebandSampleRate',fs, ...
    'Gain',txGain);

% --- RX Pluto (usb:1) ---
rx = sdrrx('Pluto', ...
    'RadioID','usb:1', ...
    'CenterFrequency',fc, ...
    'BasebandSampleRate',fs, ...
    'SamplesPerFrame',50000, ...
    'GainSource','Manual', ...
    'Gain',rxGain, ...
    'OutputDataType','double');

% Build complex tone
t = (0:N-1).'/fs;
xTone = 0.6*exp(1j*2*pi*toneHz*t);  % keep amplitude modest

% Start continuous transmit
transmitRepeat(tx, xTone);

% Capture a few frames and estimate tone power
P = zeros(1,4);
for k=1:4
    y = rx();
    P(k) = mean(abs(y).^2);
    fprintf('Frame %d power = %.3e\n', k, P(k));
end

% Stop TX and release
release(tx);
release(rx);

fprintf('Avg received power: %.3e\n', mean(P));