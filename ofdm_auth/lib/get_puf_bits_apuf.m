function R = get_puf_bits_apuf(node, C_eff, B)
% Generate B response bits from APUF by "tapping" C_eff with slight shifts.
% This avoids needing B independent challenges while staying deterministic.

m = node.m;
C_eff = logical(C_eff(:).');
if numel(C_eff) ~= m
    error('C_eff must be length m.');
end

R = false(1, B);

% Simple deterministic bit expansion: circularly rotate challenge
for i = 1:B
    C_i = circshift(C_eff, [0, i-1]);   % rotate by i-1
    R(i) = apuf_bit_from_C(C_i, node.W, node.noise_sigma);
end

R = double(R);
end
