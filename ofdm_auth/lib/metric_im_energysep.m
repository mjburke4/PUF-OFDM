function T = metric_im_energysep(Z, Atilde, N)
% Energy separation statistic over predicted active set Atilde (indices into 1..N)
e = abs(Z).^2;

mask = false(N,1);
mask(Atilde) = true;

Ea = mean(e(mask));
Ei = mean(e(~mask));
T = Ea - Ei;
end
