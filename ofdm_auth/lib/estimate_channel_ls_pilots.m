function Hhat = estimate_channel_ls_pilots(Y, Xpilot_templ, pilot_idx, N)
% LS estimate on pilots, simple linear interpolation across frequency
% Y: received FFT (N x 1)
% Xpilot_templ: expected pilot template (N x 1) (nonzero on pilots)
% pilot_idx: pilot indices
% N: FFT size

Hhat = ones(N,1);

Hp = Y(pilot_idx) ./ (Xpilot_templ(pilot_idx) + 1e-12);

% interpolate magnitude and phase separately (simple and stable)
k = pilot_idx(:);
mag = abs(Hp);
ph  = unwrap(angle(Hp));

kq = (1:N).';

mag_i = interp1(k, mag, kq, 'linear', 'extrap');
ph_i  = interp1(k, ph,  kq, 'linear', 'extrap');

Hhat = mag_i .* exp(1j*ph_i);
end
