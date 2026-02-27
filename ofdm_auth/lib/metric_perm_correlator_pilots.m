function T = metric_perm_correlator_pilots(Z, templ, pilot_idx)
% Pilot-only normalized correlator (robust scale)
% Z: equalized FFT bins (N x 1)
% templ: predicted pilot template (N x 1), nonzero on pilot bins
% pilot_idx: indices of pilot bins

zP = Z(pilot_idx);
tP = templ(pilot_idx);

den = (norm(tP)^2 + 1e-12);
T = (tP' * zP) / den;   % complex
end
