function cfg = cfg_default_dppa()
    cfg = struct();

    % Frame / auth settings
    cfg.Lsym  = 8;          % start same as your perm/IM FAR plots
    cfg.alpha = 0.75;       % frame vote acceptance ratio
    cfg.tauD  = 0.90;       % differential phase score threshold

    % PUF / challenge settings
    cfg.m = 64;
    cfg.B = 64;

    % Pilot / DPPA settings
    cfg.Npilots = 4; % for simulations %4;
    cfg.phi_choices_deg = [-25 +25]; %[-8 +8];   % start with the OTA-friendly alphabet

    % Channel / noise
    cfg.chan.type = 'awgn';

    cfg.adv_mode = 'impostor_othernode'; %'impostor_random';   % options:
    % 'impostor_random'   -> random seed / random structure
    % 'impostor_othernode'-> different enrolled node
    % 'replay_stale_nonce'-> same node, stale nonce

 %   cfgD.adv_mode = cfg.adv_mode;
end