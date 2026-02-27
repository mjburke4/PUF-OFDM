function nodes = init_node_population(U, m, varargin)
% init_node_population  Create U distinct nodes compatible with arbiter_puf_sim_node
%
% Each node is a scalar struct with field:
%   A : [m x puf_bits] matrix
%
% Optional:
%   puf_bits   : number of columns in A (latent arbiter dimensions)
%   seed       : RNG seed for reproducibility

p = inputParser;
addParameter(p,'puf_bits',128);   % internal latent dimension (can be 64,128,256)
addParameter(p,'seed',7);
parse(p,varargin{:});

puf_bits = p.Results.puf_bits;
rng(p.Results.seed);

nodes = struct([]);
for u = 1:U
    nodes(u).id = u;

    % IMPORTANT: arbiter_puf_sim_node expects this field:
    % A is challenge_len x puf_bits
    nodes(u).A = randn(m, puf_bits);

    % In init_node_population.m, inside the for u=1:U loop, add:
    nodes(u).distance_m = 5 + (200-5)*rand();   % random distance in [5,200] meters
    nodes(u).pl_d0_m = 1;                       % reference distance
    nodes(u).pl_gamma = 2.7;                    % pathloss exponent (tune)

    % Keep for reference (optional)
    nodes(u).m = m;
    nodes(u).puf_bits = puf_bits;
end
end

