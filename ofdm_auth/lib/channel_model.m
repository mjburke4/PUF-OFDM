function [y, H] = channel_model(cfg, x)
% Time-domain channel + AWGN is handled elsewhere; this returns convolution output and channel FFT.

N  = cfg.Nfft;
cp = cfg.Ncp;

if strcmpi(cfg.chan.type,'awgn')
    h = 1;
else
    Lt = cfg.chan.Ltaps;
    pow = exp(-cfg.chan.pow_exp*(0:Lt-1));
    pow = pow / sum(pow);
    h = (randn(Lt,1)+1j*randn(Lt,1))/sqrt(2);
    h = h .* sqrt(pow(:));
end

% Convolution (CP makes it circular-ish if channel shorter than CP)
y = conv(x, h);
y = y(1:numel(x)); % truncate to original length for simplicity

% Frequency response at FFT bins (assume channel shorter than CP for now)
H = fft([h; zeros(N-numel(h),1)]);
end
