function T = metric_perm_correlator(Z, templ)
% T = (templ^H Z) / ||templ||^2  (complex), caller uses |T|
den = (norm(templ)^2 + 1e-12);
T = (templ' * Z) / den;
end
