function [X, Xp, meta] = gen_ofdm_symbol(cfg)
% Generates one logical OFDM symbol X (data+pilots), and pilot-only Xp.

N = cfg.Nfft;

% Define pilot indices (simple evenly spaced over usable bins)
% For now assume fftshifted index space excluding DC
pilot_idx = round(linspace(2, N-1, cfg.Npilots));
pilot_idx = unique(pilot_idx);
pilot_idx = pilot_idx(1:cfg.Npilots);

% Known pilot values (unit magnitude)
%pilots = exp(1j * 2*pi * (0:cfg.Npilots-1)/cfg.Npilots);
pilots = exp(1j*2*pi*(0:cfg.Npilots-1)/cfg.Npilots).';


% Build X
X = zeros(N,1);
X(pilot_idx) = pilots(:);

% % Fill remaining bins with QPSK (excluding pilots)
% data_idx = setdiff(1:N, pilot_idx);
% data_bits = randi([0 1], 2*numel(data_idx), 1);
% sym = (2*data_bits(1:2:end)-1) + 1j*(2*data_bits(2:2:end)-1);
% sym = sym / sqrt(2);
% X(data_idx) = sym;

% Pilot-only vector
Xp = zeros(N,1);
Xp(pilot_idx) = X(pilot_idx);

meta.pilot_idx = pilot_idx(:);
%meta.data_idx  = data_idx(:);
meta.pilots    = pilots(:);
end
