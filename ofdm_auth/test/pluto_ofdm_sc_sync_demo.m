%% Pluto OFDM with Schmidl-Cox timing/CFO (2 Plutos)
clear; clc;

%% Settings
fc = 0.915e9;
%fc = 2.412e9;
%fc = 2.3e9;
fs = 1e6;
txGain = -2;
rxGain = 45;

Nfft = 64;
Ncp  = 16;
Nsym = 8;

Np = 16; %8;
% pilotIdx = round(linspace(2, Nfft-1, Np));
% pilotIdx = unique(pilotIdx); pilotIdx = pilotIdx(1:Np);
used = (7:58).';                 % avoid edges (tweak if needed)
pilotIdx = round(linspace(used(1), used(end), Np));
pilotIdx = unique(pilotIdx); pilotIdx = pilotIdx(1:Np);
pilots = exp(1j*2*pi*(0:Np-1)/Np).';

guard = zeros(400,1);   % bigger guard helps a lot OTA

% ----- LTF-like training symbol (frequency-domain known on used bins)
Xtr = zeros(Nfft,1);
Xtr(used) = 2*randi([0 1], numel(used), 1) - 1;   % BPSK ±1

xtr = ifft(Xtr, Nfft);
xtr_cp = [xtr(end-Ncp+1:end); xtr];

%% Build Schmidl-Cox preamble: 2 identical halves of length Nfft
% Start from random BPSK in frequency domain on even bins (classic SC)
Xsc = zeros(Nfft,1);
evenBins = 2:2:Nfft;                % even indices
bpsk = 2*randi([0 1], numel(evenBins), 1)-1;
Xsc(evenBins) = bpsk;

psc_half = ifft(Xsc, Nfft);         % time-domain length Nfft
psc = [psc_half; psc_half];         % length 2*Nfft, repeated halves
psc = 0.6 * psc / max(abs(psc)+1e-12);

%% Build OFDM payload
payload_td = zeros((Nfft+Ncp)*Nsym, 1);
for s = 1:Nsym
    X = zeros(Nfft,1);
    X(pilotIdx) = pilots;

    dataIdx = setdiff(used, pilotIdx(:));
    bits = randi([0 1], 2*numel(dataIdx), 1);
    qpsk = (2*bits(1:2:end)-1) + 1j*(2*bits(2:2:end)-1);
    qpsk = qpsk / sqrt(2);
    X(dataIdx) = qpsk;
    %dataIdx = setdiff((1:Nfft).', pilotIdx(:));
    %X(dataIdx) = 0;

    if s==1
        fprintf('TX pilot power %.3e | TX nonpilot power %.3e\n', ...
    mean(abs(X(pilotIdx)).^2), ...
    mean(abs(X(setdiff(1:Nfft,pilotIdx))).^2));
    end

    x = ifft(X, Nfft);
    x_cp = [x(end-Ncp+1:end); x];
    payload_td((s-1)*(Nfft+Ncp)+1 : s*(Nfft+Ncp)) = x_cp;
end
payload_td = 0.6 * payload_td / max(abs(payload_td)+1e-12);

%% Frame
%frame = [guard; psc; guard; payload_td; guard];
frame = [guard; psc; guard; xtr_cp; guard; payload_td; guard];
frame = 0.3 * frame / max(abs(frame)+1e-12);

%% Pluto TX/RX
tx = sdrtx('Pluto','RadioID','usb:0','CenterFrequency',fc,'BasebandSampleRate',fs,'Gain',txGain);
rx = sdrrx('Pluto','RadioID','usb:1','CenterFrequency',fc,'BasebandSampleRate',fs, ...
    'SamplesPerFrame',250000,'GainSource','Manual','Gain',rxGain,'OutputDataType','double');


%tx.transmit(frame);  % if you have the method
frameLen = numel(frame);
rx.SamplesPerFrame = frameLen + 20000;  % small padding

transmitRepeat(tx, frame);
pause(0.2);
y = rx();

release(tx);
release(rx);

%% Schmidl-Cox timing metric
L = Nfft;
ylen = numel(y);
Dmax = ylen - 2*L - 1;

P = zeros(Dmax,1);
R = zeros(Dmax,1);
for d = 1:Dmax
    a = y(d:d+L-1);
    b = y(d+L:d+2*L-1);
    %P(d) = sum(a .* conj(b));
    %R(d) = sum(abs(b).^2);
        P(d) = sum(conj(a).*b);                 % note conj(a)*b (equivalent)
    R(d) = sum(abs(a).^2) + sum(abs(b).^2); % energy of both halves
end
%M = (abs(P).^2) ./ (R.^2 + 1e-12);
M = (abs(P).^2) ./ ((0.5*R).^2 + 1e-12);    % normalize (0.5*R ~= energy of one half)
M = real(M);

fprintf('RX max amplitude %.3f\n', max(abs(y)));

% Find peak
%[~, d0] = max(M);
%fprintf('Estimated preamble start d0=%d\n', d0);

% Smooth it to reveal plateau/peak
M_s = movmean(M, 64);
% thr = 0.6 * max(M_s);              % threshold relative to peak
% minSearch = numel(guard) + 1;  % skip the first guard region
% idx = find(M_s(minSearch:end) > thr, 1, 'first');
% 
% if isempty(idx)
%     error('No Schmidl-Cox peak found. Try increasing RX gain or lowering TX gain.');
% end
% Search only in the first half of the capture (first frame)
%searchMax = min(length(M_s), 80000);  % tune if needed
%[~, d0] = max(M_s(1:searchMax));
thr = 0.6*max(M_s);
minSearch = numel(guard)+1;
idx = find(M_s(minSearch:end) > thr, 1, 'first');
d0 = idx + minSearch - 1;

%d0 = idx + minSearch - 1;
fprintf('Estimated start d0=%d\n', d0);

% CFO estimate (rad/sample)
eps_hat = angle(P(d0)) / L;
fprintf('Estimated CFO eps=%.4g rad/sample\n', eps_hat);

% CFO correct whole capture
n = (0:ylen-1).';
y_corr = y .* exp(-1j*eps_hat*n);

%% Extract first OFDM symbol after preamble+guard
preambleLen = 2*Nfft;
% Coarse location: start of training symbol (after preamble + guard)
trainStart0 = d0 + preambleLen + numel(guard);

% ---- Fine timing refinement around training start (use CP metric on training+first data symbols)
search = -Ncp:Ncp;
score = zeros(size(search));

NsCheck = 3;  % check training + 2 data symbols
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

% Now define where data starts (training symbol + guard)
dataStart = trainStart + (Nfft+Ncp) + numel(guard);

% payloadStart = d0 + preambleLen + numel(guard);
% trainStart = payloadStart;
% dataStart  = trainStart + (Nfft+Ncp) + numel(guard);  % training + guard
% 
% % bounds check
% if trainStart + (Nfft+Ncp) - 1 > numel(y_corr)
%     error('Not enough samples for training symbol. Increase SamplesPerFrame.');
% end
% 
% tr = y_corr(trainStart : trainStart + (Nfft+Ncp) - 1);
% tr_nocp = tr(Ncp+1:end);
% Ytr = fft(tr_nocp, Nfft);
% 
% % Channel estimate on used bins
% Hhat = ones(Nfft,1);
% Hhat(used) = Ytr(used) ./ (Xtr(used) + 1e-12);
% 
% % Optional smoothing across frequency
% Hhat = movmean(Hhat,5);
% 
% % --- Fine timing refinement over multiple symbols (Step 4)
% search = -Ncp:Ncp;
% score = zeros(size(search));
% 
% symIdx0 = 3;     % start checking at symbol 3 (avoid transient)
% NsCheck = 3;     % check 3 consecutive symbols
% 
% for ii = 1:numel(search)
%     s0 = payloadStart + search(ii) + (symIdx0-1)*(Nfft+Ncp);
% 
%     if s0 < 1 || (s0 + NsCheck*(Nfft+Ncp) - 1) > numel(y_corr)
%         score(ii) = -Inf;
%         continue;
%     end
% 
%     sc = 0;
%     for kk = 0:NsCheck-1
%         seg = y_corr(s0 + kk*(Nfft+Ncp) : s0 + (kk+1)*(Nfft+Ncp) - 1);
%         cp = seg(1:Ncp);
%         tail = seg(Nfft+1:Nfft+Ncp);
%         sc = sc + abs(sum(conj(cp).*tail));
%     end
%     score(ii) = sc;
% end
% 
% [~, imax] = max(score);
% shift = search(imax);
% payloadStart = payloadStart + shift;
% fprintf('Fine timing shift applied: %d samples\n', shift);
% 
% needLen = payloadStart + (Nfft+Ncp) - 1;
% if needLen > numel(y_corr)
%     error('Not enough samples after payloadStart. Increase rx.SamplesPerFrame. payloadStart=%d, needLen=%d, ylen=%d', ...
%         payloadStart, needLen, numel(y_corr));
% end


%% ---- Multi-symbol equalization and plotting (replace single-symbol block)

NsPlot = 10;      % number of OFDM symbols to accumulate
Zall = [];

for symIdx = 3:(3+NsPlot-1)

    symStart = payloadStart + (symIdx-1)*(Nfft+Ncp);

    if symStart + (Nfft+Ncp) - 1 > numel(y_corr)
        break
    end

    % local CP refine per symbol
    localSearch = -2:2;
    bestSc = -Inf; bestShift = 0;
    
    for ss = localSearch
        s_try = symStart + ss;
        if s_try < 1 || (s_try + Nfft+Ncp-1) > numel(y_corr), continue; end
        seg = y_corr(s_try : s_try + (Nfft+Ncp) - 1);
        cp = seg(1:Ncp);
        tail = seg(Nfft+1:Nfft+Ncp);
        sc = abs(sum(conj(cp).*tail));
        if sc > bestSc
            bestSc = sc;
            bestShift = ss;
        end
    end
    
    % symStart = symStart + bestShift;
    % 
    % sym = y_corr(symStart : symStart + (Nfft+Ncp) - 1);
    % sym_nocp = sym(Ncp+1:end);
    % 
    % Yk = fft(sym_nocp, Nfft);

    symStart = dataStart + (symIdx-1)*(Nfft+Ncp);
    sym = y_corr(symStart : symStart + (Nfft+Ncp) - 1);
    sym_nocp = sym(Ncp+1:end);
    
    Yk = fft(sym_nocp, Nfft);
    
    Zk = Yk ./ (Hhat + 1e-12);   % equalize with training-based estimate

    % --- Channel estimate from pilots
 
    % --- phase slope correction (SFO / residual timing)
    k = pilotIdx(:);    
    kq = (1:Nfft).';
    phi = unwrap(angle(Zk(pilotIdx).*conj(pilots)));

    pfit = polyfit(k, phi(:), 1);

    phi_corr = pfit(1)*kq + pfit(2);

    Zk = Zk .* exp(-1j*phi_corr);

    % --- common phase error correction
    phi_cpe = angle(sum(conj(pilots).*Zk(pilotIdx)));
    Zk = Zk * exp(-1j*phi_cpe);

    dataIdx = setdiff(used, pilotIdx(:));

    Zall = [Zall; Zk(dataIdx)];

end

figure;
plot(real(Zall), imag(Zall), '.');
grid on
title('Equalized QPSK constellation (multi-symbol)')
xlabel('I')
ylabel('Q')

% ---- EVM on aggregated data symbols (Zall) ----
Z = Zall(:);

% QPSK hard slicer (decide symbols first)
Shat = sign(real(Z)) + 1j*sign(imag(Z));
Shat = Shat / sqrt(2);

% Best-fit complex scalar gain (accounts for amplitude+phase offset)
g = (Shat' * Z) / (Shat' * Shat);

% EVM
evm_rms = sqrt( mean(abs(Z - g*Shat).^2) / mean(abs(g*Shat).^2) );
fprintf('EVM_rms = %.3f (%.2f dB)\n', evm_rms, 20*log10(evm_rms));

e = abs(Z - g*Shat).^2;
e_sorted = sort(e);
keep = round(0.95*numel(e));   % keep 95% best points
evm_trim = sqrt( mean(e_sorted(1:keep)) / mean(abs(g*Shat).^2) );
fprintf('EVM_trim(95%%) = %.3f (%.2f dB)\n', evm_trim, 20*log10(evm_trim));

% Remove any residual common phase rotation (best-fit)
%phi = angle(sum(Z.^2));     % for QPSK, 4th-power methods also possible
%Z = Z * exp(-1j*phi/2);

% Slice to nearest QPSK point
Zhat = sign(real(Z)) + 1j*sign(imag(Z));
Zhat = Zhat / sqrt(2);

% RMS EVM (normalized)

dataIdx = setdiff(used, pilotIdx(:));
%figure; plot(M); grid on;
%title('Schmidl-Cox timing metric M(d)'); xlabel('d'); ylabel('M(d)');

figure; plot(M_s); grid on;
title('Schmidl-Cox timing metric M(d) (smoothed)'); xlabel('d'); ylabel('M(d)');

%figure; plot(real(Y1(dataIdx)), imag(Y1(dataIdx)), '.'); grid on;
%title('Constellation after CFO correction (no EQ)'); xlabel('I'); ylabel('Q');
clipFrac = mean(abs(y) > 0.95);
fprintf('Clip fraction = %.3f\n', clipFrac);