function g = pathloss_gain(d_m, d0_m, gamma)
% pathloss_gain  Large-scale pathloss power gain (linear)
% g is power gain, so amplitude scale is sqrt(g)
g = (d0_m ./ d_m) .^ gamma;
end
