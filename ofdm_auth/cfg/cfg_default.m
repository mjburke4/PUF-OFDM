function cfg = cfg_default()
% Central config for OFDM authentication sims (edit here, everything follows)

cfg.Nfft = 64;
cfg.Ncp  = 16;

% Subcarrier indexing convention: 1..Nfft in MATLAB.
% We will treat "DC" as k = Nfft/2 + 1 for an fftshifted view if needed.
cfg.use_fftshift = true;

% Pilot tones (example: 8 pilots)
cfg.Npilots = 8;

% Modulation
cfg.mod = 'QPSK';

% Frame/voting
cfg.Lsym = 32;        % OFDM symbols per frame
cfg.alpha = 0.75;     % vote fraction threshold

% IM parameters (pilot-domain IM)
cfg.K_active_pilots = 4;  % sparsity knob: number of active pilot tones

% Permutation parameters (pilot permutation depth)
cfg.permute_pilots_only = true;
cfg.perm_depth = cfg.Npilots; % number of pilot tones being permuted (<= Npilots)

% Channel
cfg.chan.type = 'awgn';   % 'awgn' or 'rayleigh'
cfg.chan.Ltaps = 6;       % taps for rayleigh
cfg.chan.pow_exp = 0.8;   % exponential decay factor

% CSI
cfg.csi.mode = 'perfect'; % 'perfect' then later 'pilot-ls'

% Trials
cfg.MC_trials = 200;      % per SNR point
cfg.SNRdB_vec = -10:2:20;

% Thresholds (we’ll sweep in figure scripts)
%Single calibration at SNR=15 dB, L=8 (need 6/8), Vmax=3:
  %tauP=0.1800 (FAR~0.004, mean H0 votes=2.12)
  %tauIM=0.0450 (FAR~0.008, mean H0 votes=2.21)
cfg.tauP  = 0.1800;
cfg.tauIM = 0.0450;
%cfg.tauP  = 0.6;
%cfg.tauIM = 0.2;
end
