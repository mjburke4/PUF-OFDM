function T = metric_perm_correlator_activepilots(Z, templ, Atilde)
% Correlate only on predicted active pilot bins (Atilde)
zA = Z(Atilde);
tA = templ(Atilde);

den = (norm(tA)^2 + 1e-12);
T = (tA' * zA) / den;
end
