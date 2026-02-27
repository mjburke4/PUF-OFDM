function fig_frr_vs_puf_flip()

cfg = cfg_default();
cfg.chan.type = 'awgn';

% Choose operating parameters
cfg.alpha = 0.75;
L_list = [4 32];             % show short vs longer frame
SNR_list = [15];           % choose 2 representative SNRs
MC = 10000;                   % increase if you want smoother curves

% Bit flip sweep (per-bit probability)
p_list = [0 0.002 0.005 0.01 0.02 0.05];

% Nodes
U = 10; m = 64;
nodes = init_node_population(U, m, 'seed', 7);

% Fix distance so SNRref is meaningful
d0 = 10;
for u = 1:U
    nodes(u).distance_m = d0;
end
u0 = 3;

% Thresholds: choose an operating point at p_flip = 0
% (Use your existing calibration method, or set manually for now.)
% If you already have calibrated cfg.tauP/tauIM from previous work, keep them.
% Otherwise, set something reasonable to start:
if ~isfield(cfg,'tauP') || isempty(cfg.tauP);   cfg.tauP  = 0.15; end
if ~isfield(cfg,'tauIM') || isempty(cfg.tauIM); cfg.tauIM = 0.05; end

% Results arrays:
% dims: [L_idx x SNR_idx x p_idx]
FRRp = zeros(numel(L_list), numel(SNR_list), numel(p_list));
FARp = zeros(numel(L_list), numel(SNR_list), numel(p_list));
FRRi = zeros(numel(L_list), numel(SNR_list), numel(p_list));
FARi = zeros(numel(L_list), numel(SNR_list), numel(p_list));

for li = 1:numel(L_list)
    cfgL = cfg;
    cfgL.Lsym = L_list(li);
    need = ceil(cfgL.alpha * cfgL.Lsym);

    fprintf('\n==== L=%d (need %d/%d), alpha=%.2f, tauP=%.3f, tauIM=%.3f ====\n', ...
        cfgL.Lsym, need, cfgL.Lsym, cfgL.alpha, cfgL.tauP, cfgL.tauIM);

    for si = 1:numel(SNR_list)
        snr = SNR_list(si);

        for pi = 1:numel(p_list)
            p_flip = p_list(pi);

            % ---- Legit trials (H1): FRR rises with p_flip
            rejP = 0; rejI = 0;
            for t = 1:MC
                o = run_one_trial_auth(cfgL, nodes, u0, snr, true, p_flip);
                rejP = rejP + (~o.accept_perm);
                rejI = rejI + (~o.accept_im);
            end
            FRRp(li, si, pi) = rejP/MC;
            FRRi(li, si, pi) = rejI/MC;

            % ---- Impostor trials (H0): FAR (should be relatively insensitive to p_flip)
            accP = 0; accI = 0;
            for t = 1:MC
                o = run_one_trial_auth(cfgL, nodes, u0, snr, false, 0); % impostor doesn't have PUF
                accP = accP + (o.accept_perm);
                accI = accI + (o.accept_im);
            end
            FARp(li, si, pi) = accP/MC;
            FARi(li, si, pi) = accI/MC;

            fprintf('SNR= %g p=%.3f: FRR(P)=%.3f FRR(IM)=%.3f | FAR(P)=%.3f FAR(IM)=%.3f\n', ...
                snr, p_flip, FRRp(li,si,pi), FRRi(li,si,pi), FARp(li,si,pi), FARi(li,si,pi));
        end
    end
end

% ---- Plot FRR vs p_flip (Permutation)
figure; hold on; grid on; set(gca,'YScale','log');
for li = 1:numel(L_list)
    for si = 1:numel(SNR_list)
        plot(p_list, max(squeeze(FRRp(li,si,:))', 1/MC), '-o', 'LineWidth', 2);
    end
end
xlabel('PUF bit-flip probability p_{flip}');
ylabel('FRR (log scale)');
title(sprintf('Permutation: FRR vs PUF Bit-Flip (\\alpha=%.2f, MC=%d)', cfg.alpha, MC));

leg = strings(1, numel(L_list)*numel(SNR_list));
k = 1;
for li = 1:numel(L_list)
    for si = 1:numel(SNR_list)
        leg(k) = sprintf('L=%d, SNR= %g dB', L_list(li), SNR_list(si));
        k = k + 1;
    end
end
legend(leg, 'Location','best');

% ---- Plot FRR vs p_flip (IM)
figure; hold on; grid on; set(gca,'YScale','log');
for li = 1:numel(L_list)
    for si = 1:numel(SNR_list)
        plot(p_list, max(squeeze(FRRi(li,si,:))', 1/MC), '-o', 'LineWidth', 2);
    end
end
xlabel('PUF bit-flip probability p_{flip}');
ylabel('FRR (log scale)');
title(sprintf('IM: FRR vs PUF Bit-Flip (\\alpha=%.2f, MC=%d)', cfg.alpha, MC));
legend(leg, 'Location','best');

end