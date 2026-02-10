function Hhat = estimate_channel(cfg, Htrue, Zpilot, Xpilot, pilot_idx)
% Perfect first; later add LS estimate.
if strcmpi(cfg.csi.mode,'perfect')
    Hhat = Htrue;
else
    Hhat = ones(size(Htrue));
    Hhat(pilot_idx) = Zpilot(pilot_idx) ./ Xpilot(pilot_idx);
    % crude interpolation left for later
end
end
