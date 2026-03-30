function out = run_one_trial_dppa(cfg, nodes, u0, SNRref_dB, is_legit, p_flip)
% run_one_trial_dppa
% One frame of DPPA authentication using a differential pilot-phase detector.
%
% Inputs:
%   cfg         : config struct from cfg_default_dppa()
%   nodes       : node population from init_node_population(...)
%   u0          : claimed legitimate node index
%   SNRref_dB   : reference SNR in dB
%   is_legit    : true = legitimate transmitter, false = impostor
%   p_flip      : optional PUF bit flip probability for TX side
%
% Output:
%   out contains per-symbol scores, votes, and frame decision

if nargin < 6, p_flip = 0; end

L = cfg.Lsym;
B = cfg.B;
m = cfg.m;
Npilots = cfg.Npilots;

% Frame-level public values
C_seed = randi([0 1], 1, m);
nonce_frame = randi([0 1], 1, 32);

% Base pilot vector (simple pilot-domain abstraction)
p_base = ones(Npilots,1);

% Choose TX identity
% if is_legit
%     node_tx = nodes(u0);
% else
%     % Impostor: choose a different node if available, else random phase behavior
%     if numel(nodes) >= 2
%         imp_choices = setdiff(1:numel(nodes), u0);
%         u_imp = imp_choices(randi(numel(imp_choices)));
%         node_tx = nodes(u_imp);
%     else
%         node_tx = nodes(u0); % fallback if only one node exists
%     end
% end

u_tx = u0;

if ~is_legit
    switch cfg.adv_mode
        case 'impostor_othernode'
            cand = setdiff(1:numel(nodes), u0);
            u_tx = cand(randi(numel(cand)));
        case {'impostor_random','replay_stale_nonce'}
            u_tx = u0;
        otherwise
            error('Unknown cfg.adv_mode: %s', cfg.adv_mode);
    end
end

node_tx = nodes(u_tx);
% Generate received pilot observations per symbol
u_meas = zeros(1,L);
phi_tx_deg = zeros(1,L);
phiH1_deg  = zeros(1,L);
phiH0_deg  = zeros(1,L);

snr_lin = 10^(SNRref_dB/10);
sigma2 = 1/snr_lin;  % simple pilot-domain AWGN normalization


for ell = 1:L
    % ----- TX effective challenge -----
    nonce_bits_tx = nonce_with_symbol(nonce_frame, ell, 16);
    C_eff_tx      = mix_challenge_nonce(C_seed, nonce_bits_tx, m);

    if is_legit
        R_tx = get_puf_bits(nodes(u0), C_eff_tx, B);

    else
        switch cfg.adv_mode
            case 'impostor_random'
                % random phase seed impostor
                seed_tx = uint32(randi([0, 2^32-1]));
                phi_tx_deg(ell) = phase_from_seed_dppa(seed_tx, cfg.phi_choices_deg);
                phi_tx = deg2rad(phi_tx_deg(ell));
                p_tx = p_base * exp(1j*phi_tx);

                w = sqrt(sigma2/2) * (randn(Npilots,1) + 1j*randn(Npilots,1));
                z = p_tx + w;
                u_meas(ell) = (p_base' * z) / (norm(p_base)^2 + 1e-12);

                % H1 still assumes claimed node u0 with correct nonce
                nonce_bits_H1 = nonce_with_symbol(nonce_frame, ell, 16);
                C_eff_H1      = mix_challenge_nonce(C_seed, nonce_bits_H1, m);
                R_H1          = get_puf_bits(nodes(u0), C_eff_H1, B);
                seed_H1       = hash_bits_u32(R_H1);
                phiH1_deg(ell)= phase_from_seed_dppa(seed_H1, cfg.phi_choices_deg);

                % Optional wrong-nonce reference for separate replay analysis
                nonce_frame_H0 = [1 - nonce_frame(1), nonce_frame(2:end)];
                nonce_bits_H0  = nonce_with_symbol(nonce_frame_H0, ell, 16);
                C_eff_H0       = mix_challenge_nonce(C_seed, nonce_bits_H0, m);
                R_H0           = get_puf_bits(nodes(u0), C_eff_H0, B);
                seed_H0        = hash_bits_u32(R_H0);
                phiH0_deg(ell) = phase_from_seed_dppa(seed_H0, cfg.phi_choices_deg);
                continue

            case 'impostor_othernode'
                R_tx = get_puf_bits(node_tx, C_eff_tx, B);

            case 'replay_stale_nonce'
                nonce_frame_tx = [1 - nonce_frame(1), nonce_frame(2:end)];
                nonce_bits_tx  = nonce_with_symbol(nonce_frame_tx, ell, 16);
                C_eff_tx       = mix_challenge_nonce(C_seed, nonce_bits_tx, m);
                R_tx           = get_puf_bits(nodes(u0), C_eff_tx, B);

            otherwise
                error('Unknown cfg.adv_mode: %s', cfg.adv_mode);
        end
    end

    if p_flip > 0
        flipmask = rand(size(R_tx)) < p_flip;
        R_tx = xor(logical(R_tx), flipmask);
        R_tx = double(R_tx);
    end

    seed_tx = hash_bits_u32(R_tx);
    phi_tx_deg(ell) = phase_from_seed_dppa(seed_tx, cfg.phi_choices_deg);
    phi_tx = deg2rad(phi_tx_deg(ell));

    p_tx = p_base * exp(1j*phi_tx);

    % AWGN observation
    w = sqrt(sigma2/2) * (randn(Npilots,1) + 1j*randn(Npilots,1));
    z = p_tx + w;

    u_meas(ell) = (p_base' * z) / (norm(p_base)^2 + 1e-12);

    % ----- H1: verifier assumes claimed node u0 with correct nonce -----
    nonce_bits_H1 = nonce_with_symbol(nonce_frame, ell, 16);
    C_eff_H1      = mix_challenge_nonce(C_seed, nonce_bits_H1, m);
    R_H1          = get_puf_bits(nodes(u0), C_eff_H1, B);
    seed_H1       = hash_bits_u32(R_H1);
    phiH1_deg(ell)= phase_from_seed_dppa(seed_H1, cfg.phi_choices_deg);

    % ----- Optional replay / stale-nonce reference -----
    nonce_frame_H0 = [1 - nonce_frame(1), nonce_frame(2:end)];
    nonce_bits_H0  = nonce_with_symbol(nonce_frame_H0, ell, 16);
    C_eff_H0       = mix_challenge_nonce(C_seed, nonce_bits_H0, m);
    R_H0           = get_puf_bits(nodes(u0), C_eff_H0, B);
    seed_H0        = hash_bits_u32(R_H0);
    phiH0_deg(ell) = phase_from_seed_dppa(seed_H0, cfg.phi_choices_deg);
end

% Differential detector
Td_H1 = zeros(1, L-1);
Td_H0 = zeros(1, L-1);
dErr_H1_deg = zeros(1, L-1);
dErr_H0_deg = zeros(1, L-1);

for ell = 2:L
    d_meas = u_meas(ell) * conj(u_meas(ell-1));
    dphi_meas = angle(d_meas);

    dphi_H1 = deg2rad(phiH1_deg(ell) - phiH1_deg(ell-1));
    dphi_H0 = deg2rad(phiH0_deg(ell) - phiH0_deg(ell-1));

    e1 = angle(exp(1j*(dphi_meas - dphi_H1)));
    e0 = angle(exp(1j*(dphi_meas - dphi_H0)));

    Td_H1(ell-1) = cos(e1);
    Td_H0(ell-1) = cos(e0);

    dErr_H1_deg(ell-1) = rad2deg(e1);
    dErr_H0_deg(ell-1) = rad2deg(e0);
end

% Vote-based decisions
Ldiff = numel(Td_H1);
Vd_H1 = sum(Td_H1 >= cfg.tauD);
Vd_H0 = sum(Td_H0 >= cfg.tauD);

out.Vd_H1 = Vd_H1;
out.Vd_H0 = Vd_H0;
out.accept_dppa = (Vd_H1 >= ceil(cfg.alpha * Ldiff));   % for legit/impostor trial bookkeeping
out.accept_h0   = (Vd_H0 >= ceil(cfg.alpha * Ldiff));   % wrong-hypothesis acceptance
out.Td_H1 = Td_H1;
out.Td_H0 = Td_H0;
out.dErr_H1_deg = dErr_H1_deg;
out.dErr_H0_deg = dErr_H0_deg;
out.phi_tx_deg = phi_tx_deg;
out.phiH1_deg = phiH1_deg;
out.phiH0_deg = phiH0_deg;
end

