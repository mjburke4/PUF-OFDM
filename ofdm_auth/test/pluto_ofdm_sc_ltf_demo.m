%% pluto_ofdm_sc_ltf_demo.m
% Two-Pluto OTA OFDM bring-up with:
% - Schmidl-Cox timing/CFO
% - LTF-like training symbol (dense channel estimate)
% - Fine timing around training
% - Pilot CPE correction
% - Multi-symbol constellation + EVM
%
% TX: Pluto usb:0
% RX: Pluto usb:1

clear; close all;

%% ---------------- User settings ----------------
fc = 915e6;          % center frequency (use 915 MHz for clean band)
fc = 2.437e9;
fs = 1e6;            % sample rate
txGain = -10;        % adjust as needed
rxGain = 45;         % adjust as needed (avoid saturation)

Nfft = 64;
Ncp  = 16;

% Use subcarriers (avoid DC/edges)
used = (7:58).';      % 52 bins like 802.11a-style used region

% Pilots (more pilots = better HW equalization)
Np = 16;
pilotIdx = round(linspace(used(1), used(end), Np));
pilotIdx = unique(pilotIdx); pilotIdx = pilotIdx(1:Np);
pilots = exp(1j*2*pi*(0:Np-1)/Np).';    % fixed distinct pilots

% Data indices are used minus pilots
dataIdx = setdiff(used, pilotIdx(:));

NsymPayload = 20;   % payload OFDM symbols
NsPlot = 12;        % number of symbols to accumulate for constellation
symIdxStart = 4;    % start plotting from symbol #4 (avoid transient)

guardLen = 400;     % samples
guard = zeros(guardLen,1);

% Capture length
% We'll compute frameLen after building frame and set rx accordingly.
extraPad = 20000;

% Fine timing refinement around training symbol
NsCheck = 3;         % training + 2 data symbols
% ------------------------------------------------


%% ---------------- Build Schmidl-Cox preamble ----------------
% Classic Schmidl-Cox: populate only even bins with BPSK -> time domain has two equal halves.
Xsc = zeros(Nfft,1);
evenBins = 2:2:Nfft;
bpsk_sc = 2*randi([0 1], numel(evenBins), 1)-1;
Xsc(evenBins) = bpsk_sc;

psc_half = ifft(Xsc, Nfft);
psc = [psc_half; psc_half];  % length 2*Nfft
psc = 0.5 * psc / max(abs(psc)+1e-12);

%% ---------------- Build LTF-like training symbol ----------------
Xtr = zeros(Nfft,1);
Xtr(used) = 2*randi([0 1], numel(used), 1) - 1;  % BPSK ±1 on used bins
xtr = ifft(Xtr, Nfft);
xtr_cp = [xtr(end-Ncp+1:end); xtr];
xtr_cp = 0.5 * xtr_cp / max(abs(xtr_cp)+1e-12);

%% ---------------- Build OFDM payload ----------------
payload_td = zeros((Nfft+Ncp)*NsymPayload, 1);

for s = 1:NsymPayload
    X = zeros(Nfft,1);

    % Pilots
    X(pilotIdx) = pilots;

    % QPSK data on used bins (excluding pilots)
    bits = randi([0 1], 2*numel(dataIdx), 1);
    qpsk = (2*bits(1:2:end)-1) + 1j*(2*bits(2:2:end)-1);
    qpsk = qpsk / sqrt(2);
    X(dataIdx) = qpsk;

    if s == 1
        fprintf('TX pilot power %.3e | TX nonpilot(used,data) power %.3e\n', ...
            mean(abs(X(pilotIdx)).^2), mean(abs(X(dataIdx)).^2));
    end

    x = ifft(X, Nfft);
    x_cp = [x(end-Ncp+1:end); x];

    payload_td((s-1)*(Nfft+Ncp)+1 : s*(Nfft+Ncp)) = x_cp;
end
payload_td = 0.5 * payload_td / max(abs(payload_td)+1e-12);

%% ---------------- Assemble frame ----------------
% Frame: guard | SC preamble | guard | training | guard | payload | guard
frame = [guard; psc; guard; xtr_cp; guard; payload_td; guard];
frame = 0.3 * frame / max(abs(frame)+1e-12);   % keep conservative amplitude
frameLen = numel(frame);

fprintf('Frame length = %d samples\n', frameLen);

%% ---------------- Pluto TX/RX ----------------
tx = sdrtx('Pluto', ...
    'RadioID','usb:0', ...
    'CenterFrequency',fc, ...
    'BasebandSampleRate',fs, ...
    'Gain',txGain);

rx = sdrrx('Pluto', ...
    'RadioID','usb:1', ...
    'CenterFrequency',fc, ...
    'BasebandSampleRate',fs, ...
    'SamplesPerFrame', frameLen + extraPad, ...
    'GainSource','Manual', ...
    'Gain',rxGain, ...
    'OutputDataType','double');

%% ---------------- Transmit and capture ----------------
transmitRepeat(tx, frame);
pause(0.2);
y = rx();

release(tx);
release(rx);

fprintf('RX max amplitude %.3f\n', max(abs(y)));
clipFrac = mean(abs(y) > 0.95);
fprintf('Clip fraction = %.3f\n', clipFrac);

%% ---------------- Schmidl-Cox timing/CFO ----------------
L = Nfft;
ylen = numel(y);
Dmax = ylen - 2*L - 1;

P = zeros(Dmax,1);
R = zeros(Dmax,1);

for d = 1:Dmax
    a = y(d:d+L-1);
    b = y(d+L:d+2*L-1);
    P(d) = sum(conj(a).*b);
    R(d) = sum(abs(a).^2) + sum(abs(b).^2);
end

M = (abs(P).^2) ./ ((0.5*R).^2 + 1e-12);
M = real(M);
M_s = movmean(M, 64);

% Pick first strong region after initial guard
thr = 0.6 * max(M_s);
minSearch = guardLen + 1;
idx = find(M_s(minSearch:end) > thr, 1, 'first');
if isempty(idx)
    error('No Schmidl-Cox detection. Increase RX gain or TX gain.');
end
d0 = idx + minSearch - 1;
fprintf('Estimated start d0=%d\n', d0);

eps_hat = angle(P(d0)) / L;
fprintf('Estimated CFO eps=%g rad/sample\n', eps_hat);

% CFO correct
n = (0:ylen-1).';
y_corr = y .* exp(-1j*eps_hat*n);

figure; plot(M_s); grid on;
title('Schmidl-Cox timing metric M(d) (smoothed)');
xlabel('d'); ylabel('M(d)');

%% ---------------- Coarse training start ----------------
preambleLen = 2*Nfft;
trainStart0 = d0 + preambleLen + guardLen;

%% ---------------- Fine timing around training start (multi-symbol CP metric) ----------------
search = -Ncp:Ncp;
score = zeros(size(search));

for ii = 1:numel(search)
    s0 = trainStart0 + search(ii);
    if s0 < 1 || (s0 + NsCheck*(Nfft+Ncp) - 1) > numel(y_corr)
        score(ii) = -Inf; continue;
    end

    sc = 0;
    for kk = 0:NsCheck-1
        seg = y_corr(s0 + kk*(Nfft+Ncp) : s0 + (kk+1)*(Nfft+Ncp) - 1);
        cp = seg(1:Ncp);
        tail = seg(Nfft+1:Nfft+Ncp);
        sc = sc + abs(sum(conj(cp).*tail));
    end
    score(ii) = sc;
end

[~, imax] = max(score);
shift = search(imax);
trainStart = trainStart0 + shift;
fprintf('Fine timing shift applied: %d samples\n', shift);

% Data starts after: training symbol + guard
dataStart = trainStart + (Nfft+Ncp) + guardLen;

%% ---------------- Extract training symbol and estimate channel ----------------
if trainStart + (Nfft+Ncp) - 1 > numel(y_corr)
    error('Not enough samples for training symbol. Increase SamplesPerFrame.');
end

tr = y_corr(trainStart : trainStart + (Nfft+Ncp) - 1);
tr_nocp = tr(Ncp+1:end);
Ytr = fft(tr_nocp, Nfft);

Hhat_tr = ones(Nfft,1);
Hhat_tr(used) = Ytr(used) ./ (Xtr(used) + 1e-12);
Hhat_tr = movmean(Hhat_tr, 5);

fprintf('Mean |Hhat_tr(used)| = %.3f\n', mean(abs(Hhat_tr(used))));

%% ---------------- Demod multiple payload symbols ----------------
Zall = [];
evm_list = [];

% local CP refine per symbol
localSearch = -2:2;

for symIdx = symIdxStart:(symIdxStart+NsPlot-1)

    symStart = dataStart + (symIdx-1)*(Nfft+Ncp);

    if symStart + (Nfft+Ncp) - 1 > numel(y_corr)
        break;
    end

    % Local CP refine
    bestSc = -Inf; bestShift = 0;
    for ss = localSearch
        s_try = symStart + ss;
        if s_try < 1 || (s_try + (Nfft+Ncp) - 1) > numel(y_corr), continue; end
        seg = y_corr(s_try : s_try + (Nfft+Ncp) - 1);
        cp = seg(1:Ncp);
        tail = seg(Nfft+1:Nfft+Ncp);
        sc = abs(sum(conj(cp).*tail));
        if sc > bestSc
            bestSc = sc;
            bestShift = ss;
        end
    end
    symStart = symStart + bestShift;

    sym = y_corr(symStart : symStart + (Nfft+Ncp) - 1);
    sym_nocp = sym(Ncp+1:end);
    Yk = fft(sym_nocp, Nfft);

    % Equalize with training-based channel estimate
    Zk = Yk ./ (Hhat_tr + 1e-12);

    % CPE correction using pilots
    phi_cpe = angle(sum(conj(pilots) .* Zk(pilotIdx)));
    Zk = Zk * exp(-1j*phi_cpe);

    % Collect data bins only
    Zsym = Zk(dataIdx);
    Zall = [Zall; Zsym];

    % Per-symbol EVM estimate (optional)
    Shat = sign(real(Zsym)) + 1j*sign(imag(Zsym));
    Shat = Shat / sqrt(2);
    g = (Shat' * Zsym) / (Shat' * Shat);
    evm_sym = sqrt(mean(abs(Zsym - g*Shat).^2) / mean(abs(g*Shat).^2));
    evm_list(end+1) = evm_sym; %#ok<AGROW>
end

%% ---------------- Plot constellation ----------------
figure;
plot(real(Zall), imag(Zall), '.'); grid on;
title('Equalized QPSK constellation (multi-symbol)');
xlabel('I'); ylabel('Q');

%% ---------------- EVM on aggregated symbols ----------------
Z = Zall(:);

Shat = sign(real(Z)) + 1j*sign(imag(Z));
Shat = Shat / sqrt(2);

g = (Shat' * Z) / (Shat' * Shat);

evm_rms = sqrt(mean(abs(Z - g*Shat).^2) / mean(abs(g*Shat).^2));
fprintf('EVM_rms = %.3f (%.2f dB)\n', evm_rms, 20*log10(evm_rms));

e = abs(Z - g*Shat).^2;
e_sorted = sort(e);
keep = round(0.95*numel(e));
evm_trim = sqrt(mean(e_sorted(1:keep)) / mean(abs(g*Shat).^2));
fprintf('EVM_trim(95%%) = %.3f (%.2f dB)\n', evm_trim, 20*log10(evm_trim));

fprintf('EVM per symbol: mean=%.3f, min=%.3f, max=%.3f\n', ...
    mean(evm_list), min(evm_list), max(evm_list));