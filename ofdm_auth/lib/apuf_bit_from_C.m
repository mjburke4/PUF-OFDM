function r = apuf_bit_from_C(C_bits, W, noise_sigma)
% Raw Arbiter PUF bit: r = sign(Phi*W + n) >= 0
% C_bits: 0/1 row vector length m
% W: (m+1)x1
% noise_sigma: std dev of additive Gaussian noise on decision variable

if nargin < 3, noise_sigma = 0; end

b = 1 - 2*double(C_bits(:).');     % 0->+1, 1->-1
p = fliplr(cumprod(fliplr(b),2));  % suffix product
Phi = [1, p];                      % 1 x (m+1)

d = Phi * W;
if noise_sigma > 0
    d = d + noise_sigma * randn();
end
r = (d >= 0);
end
