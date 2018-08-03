addpath '.\Wang Li code - original download';

% Example: Call MCMC sampler for multiple graphs on simulated data with 3
% subgroups on 20 variables

% p is number of variables
p = 20;

% K is number of subgroups
K = 3;

% Set up a precision matrix for each group. Here they are all the same.
Omega1 = toeplitz([1, 0.5, 0.4, zeros(1, p - 3)]);
Omega2 = Omega1;
Omega3 = Omega1;

% True covariance matrix is inverse of precision matrix
Cov1_True = inv(Omega1);
Cov2_True = inv(Omega2);
Cov3_True = inv(Omega3);

% Simulate data with sample size 100 per group
n = 100;

% Random normal data with mean 0 and given covariance and sample size
X1 = rMNorm(zeros(p, 1), Cov1_True, n)';
X2 = rMNorm(zeros(p, 1), Cov2_True, n)';
X3 = rMNorm(zeros(p, 1), Cov3_True, n)';

% X'X matrix
S1 = X1' * X1;
S2 = X2' * X2;
S3 = X3' * X3;

% Number of MCMC iterations before burnin
burnin  = 10000;

% Number of MCMC iterations after burnin
nmc = 20000;

%% Setting Hyperparameters %%
h = 50^2;
v0=0.02^2;
v1=h*v0;
lambda = 1;
pii = 2/(p-1);
pii = [pii pii pii pii];
% Graph Specific Hyperparameters %
V0_1 = v0*ones(p);
V1_1 = v1*ones(p);
V0 = cat(3,V0_1, V0_1,V0_1, V0_1);
V1 = cat(3, V1_1, V1_1,V1_1, V1_1);

% Prior parameters for gamma slab of mixture prior
alpha = 4;
beta = 5;

% Parameters for prior on nu which affects graph sparsity
a = 1;
b = 4;

% Parameter for Bernoulli prior on indicator of graph relatedness
w = .9;

% Initial value for precision matrix for each group
InitialC = eye(p);

% Initial values for Theta and nu
Theta = zeros(K);
nu = zeros(p, p) - 1;

C = cat(3, InitialC, InitialC, InitialC);
InitialSig = inv(InitialC);
Sig = cat(3, InitialSig, InitialSig, InitialSig);

% Call MCMC sampler
[C_save, Sig_save, adj_save, Theta_save, ar_gamma, ar_theta, nu_save, ar_nu] = ...
    MCMC_multiple_graphs_SSVS_Final(Theta, Sig, V0, V1, lambda, pii, repmat(n, 1, K), cat(3,S1, S2, S3), ...
    C, nu, alpha, beta, a, b, w, burnin, nmc, true);

% PPIs for Theta (graph similarity measure)
ppi_theta = mean(Theta_save ~= 0, 3);

% Edge PPIs for each graph
ppi_edges = mean(adj_save, 4);

% Get 95% credible intervals for omega (precision matrix)
CI_omega_lower = quantile(C_save, 0.025, 4);
CI_omega_upper = quantile(C_save, 0.975, 4);