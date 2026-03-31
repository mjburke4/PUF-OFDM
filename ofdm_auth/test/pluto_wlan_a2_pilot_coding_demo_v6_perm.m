clear; clc; %close all;

Ntrials = 20;
fc_list = [350e6];
tau_p   = 0.95;    % starter value, tune later
alpha   = 0.6; %0.75;

results = repmat(struct( ...
    'trial', [], ...
    'fc', [], ...
    'ber', [], ...
    'numErr', [], ...
    'Lbits', [], ...
    'TpH1_mean', [], ...
    'TpH0_mean', [], ...
    'VH1', [], ...
    'VH0', [], ...
    'Lsym', [], ...
    'Lauth', [], ...
    'acceptH1', [], ...
    'acceptH0', [], ...
    'pktStart', [], ...
    'coarseCFO', [], ...
    'fineCFO', [], ...
    'rxMaxAmp', [], ...
    'Tp_H1', [], ...
    'Tp_H0', [], ...
    'ok', false, ...
    'errmsg', ''), Ntrials, 1);

fprintf('\n===== STARTING CONSTRAINED PERM OTA RUN =====\n');

for k = 1:Ntrials
    fc_this = fc_list(mod(k-1, numel(fc_list)) + 1);

    fprintf('\n----------------------------------------\n');
    fprintf('Trial %d / %d   fc = %.3f MHz\n', k, Ntrials, fc_this/1e6);
    fprintf('----------------------------------------\n');

    try
        res = run_one_perm_trial(fc_this, tau_p, alpha);

        % Test code section changing Lsym
        %res.Lsym = 64;

        results(k).trial      = k;
        results(k).fc         = res.fc;
        results(k).ber        = res.ber;
        results(k).numErr     = res.numErr;
        results(k).Lbits      = res.Lbits;
        results(k).TpH1_mean  = res.TpH1_mean;
        results(k).TpH0_mean  = res.TpH0_mean;
        results(k).VH1        = res.VH1;
        results(k).VH0        = res.VH0;
        results(k).Lsym       = res.Lsym;
        results(k).Lauth      = res.Lauth;
        results(k).acceptH1   = res.acceptH1;
        results(k).acceptH0   = res.acceptH0;
        results(k).pktStart   = res.pktStart;
        results(k).coarseCFO  = res.coarseCFO;
        results(k).fineCFO    = res.fineCFO;
        results(k).rxMaxAmp   = res.rxMaxAmp;
        results(k).Tp_H1      = res.Tp_H1;
        results(k).Tp_H0      = res.Tp_H0;
        results(k).ok         = true;
        results(k).errmsg     = '';

        fprintf('Trial %d summary: BER=%.3e, Tp_H1=%.3f, Tp_H0=%.3f, acc(H1/H0)=%d/%d\n', ...
            k, res.ber, res.TpH1_mean, res.TpH0_mean, res.acceptH1, res.acceptH0);

    catch ME
        results(k).trial = k;
        results(k).fc    = fc_this;
        results(k).ok    = false;
        results(k).errmsg = ME.message;

        fprintf(2, 'Trial %d FAILED: %s\n', k, ME.message);
    end
end

fprintf('\n===== MULTI-TRIAL RUN COMPLETE =====\n');

ok_idx = find([results.ok]);

if isempty(ok_idx)
    error('No successful trials completed.');
end

ber_vec   = [results(ok_idx).ber];
TpH1_vec  = [results(ok_idx).TpH1_mean];
TpH0_vec  = [results(ok_idx).TpH0_mean];
accH1_vec = [results(ok_idx).acceptH1];
accH0_vec = [results(ok_idx).acceptH0];
Lsym_vec  = [results(ok_idx).Lsym];
Lauth_vec = [results(ok_idx).Lauth];

fprintf('\n===== SUMMARY OVER %d SUCCESSFUL TRIALS =====\n', numel(ok_idx));
fprintf('Mean BER               = %.3e\n', mean(ber_vec));
fprintf('Median BER             = %.3e\n', median(ber_vec));
fprintf('Mean Tp_H1             = %.3f\n', mean(TpH1_vec));
fprintf('Mean Tp_H0             = %.3f\n', mean(TpH0_vec));
fprintf('H1 acceptance rate     = %.3f\n', mean(accH1_vec));
fprintf('H0 acceptance rate     = %.3f\n', mean(accH0_vec));
fprintf('Mean payload Lsym      = %.1f\n', mean(Lsym_vec));
fprintf('Mean auth-window Lauth = %.1f\n', mean(Lauth_vec));