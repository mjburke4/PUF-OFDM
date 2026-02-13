function nodes = init_node_population_apuf(U, m, varargin)
% nodes = init_node_population_apuf(U, m, 'noise_sigma', 0.0, 'seed', 1)
% Creates U nodes, each with distinct Arbiter PUF weight vector W (size m+1).

p = inputParser;
addParameter(p,'noise_sigma',0.0);     % default: perfect PUF
addParameter(p,'seed',1);
parse(p,varargin{:});

rng(p.Results.seed);

nodes = struct([]);
for u = 1:U
    nodes(u).id = u;
    nodes(u).m  = m;
    % Typical APUF model: random weights (Gaussian)
    nodes(u).W = randn(m+1,1);   % distinct per node
    nodes(u).noise_sigma = p.Results.noise_sigma;
end
end
